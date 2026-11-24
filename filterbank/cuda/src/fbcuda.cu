#include <cuda_runtime.h>
#include <cufft.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>

#include "fbcuda.h"

static thread_local std::string g_err;

#define CK(x)                                                                 \
    do {                                                                      \
        cudaError_t e_ = (x);                                                 \
        if (e_ != cudaSuccess) {                                              \
            g_err = std::string(#x) + ": " + cudaGetErrorString(e_);          \
            return -1;                                                        \
        }                                                                     \
    } while (0)

template <typename T>
struct DevBuf {
    T *p = nullptr;
    ~DevBuf() {
        if (p) cudaFree(p);
    }
    cudaError_t alloc(size_t n) { return cudaMalloc(&p, n * sizeof(T)); }
};

__global__ void transpose_kernel(const float *in, float *out, int nsamp, int nchan) {
    __shared__ float tile[32][33];
    int c = blockIdx.x * 32 + threadIdx.x;
    int t = blockIdx.y * 32 + threadIdx.y;
    for (int j = 0; j < 32; j += 8) {
        if (c < nchan && t + j < nsamp) tile[threadIdx.y + j][threadIdx.x] = in[(size_t)(t + j) * nchan + c];
    }
    __syncthreads();
    c = blockIdx.x * 32 + threadIdx.y;
    t = blockIdx.y * 32 + threadIdx.x;
    for (int j = 0; j < 32; j += 8) {
        if (c + j < nchan && t < nsamp) out[(size_t)(c + j) * nsamp + t] = tile[threadIdx.x][threadIdx.y + j];
    }
}

__global__ void dedisperse_kernel(const float *tr, int nsamp, int nchan, const int *delays,
                                  const unsigned char *kill, int nout, float *out) {
    extern __shared__ int sdel[];
    int d = blockIdx.y;
    for (int c = threadIdx.x; c < nchan; c += blockDim.x) {
        int v = delays[(size_t)d * nchan + c];
        sdel[c] = (kill && kill[c]) ? -1 : v;
    }
    __syncthreads();
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= nout) return;
    float acc = 0.0f;
    for (int c = 0; c < nchan; ++c) {
        int off = sdel[c];
        if (off >= 0) acc += tr[(size_t)c * nsamp + t + off];
    }
    out[(size_t)d * nout + t] = acc;
}

__global__ void zero_dm_kernel(float *data, int nsamp, int nchan, const float *chan_mean) {
    int t = blockIdx.x;
    if (t >= nsamp) return;
    __shared__ float part[256];
    float s = 0.0f;
    for (int c = threadIdx.x; c < nchan; c += blockDim.x) s += data[(size_t)t * nchan + c];
    part[threadIdx.x] = s;
    __syncthreads();
    for (int k = blockDim.x / 2; k > 0; k >>= 1) {
        if (threadIdx.x < k) part[threadIdx.x] += part[threadIdx.x + k];
        __syncthreads();
    }
    float m = part[0] / nchan;
    for (int c = threadIdx.x; c < nchan; c += blockDim.x) {
        data[(size_t)t * nchan + c] += (chan_mean ? chan_mean[c] : 0.0f) - m;
    }
}

__global__ void channel_mean_kernel(const float *data, int nsamp, int nchan, float *mean) {
    int c = blockIdx.x * blockDim.x + threadIdx.x;
    if (c >= nchan) return;
    double s = 0.0;
    for (int t = 0; t < nsamp; ++t) s += data[(size_t)t * nchan + c];
    mean[c] = (float)(s / nsamp);
}

__global__ void boxcar_kernel(const double *cs, int ndm, int nout, const int *widths, int nwidth, const double *noise,
                              float *best_snr, int *best_width) {
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    int d = blockIdx.y;
    if (t >= nout) return;
    const double *row = cs + (size_t)d * (nout + 1);
    float best = -1e30f;
    int bw = 0;
    for (int k = 0; k < nwidth; ++k) {
        int w = widths[k];
        if (t + w > nout) break;
        const double *nz = noise + ((size_t)d * nwidth + k) * 2;
        float s = (float)((row[t + w] - row[t] - nz[0]) * nz[1]);
        if (s > best) {
            best = s;
            bw = w;
        }
    }
    best_snr[(size_t)d * nout + t] = best;
    best_width[(size_t)d * nout + t] = bw;
}

extern "C" int fbcuda_device_count(void) {
    int n = 0;
    if (cudaGetDeviceCount(&n) != cudaSuccess) return 0;
    return n;
}

extern "C" const char *fbcuda_last_error(void) { return g_err.c_str(); }

extern "C" int fbcuda_set_device(int device) {
    CK(cudaSetDevice(device));
    return 0;
}

extern "C" int fbcuda_dedisperse(const float *data, int nsamp, int nchan, const int *delays, int ndm,
                                 const unsigned char *killmask, int nout, float *out) {
    DevBuf<float> d_in, d_tr, d_out;
    DevBuf<int> d_del;
    DevBuf<unsigned char> d_kill;
    size_t n = (size_t)nsamp * nchan;
    CK(d_in.alloc(n));
    CK(d_tr.alloc(n));
    CK(d_out.alloc((size_t)ndm * nout));
    CK(d_del.alloc((size_t)ndm * nchan));
    CK(cudaMemcpy(d_in.p, data, n * sizeof(float), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_del.p, delays, (size_t)ndm * nchan * sizeof(int), cudaMemcpyHostToDevice));
    if (killmask) {
        CK(d_kill.alloc(nchan));
        CK(cudaMemcpy(d_kill.p, killmask, nchan, cudaMemcpyHostToDevice));
    }
    dim3 tb(32, 8), tg((nchan + 31) / 32, (nsamp + 31) / 32);
    transpose_kernel<<<tg, tb>>>(d_in.p, d_tr.p, nsamp, nchan);
    CK(cudaGetLastError());
    int threads = 256;
    dim3 g((nout + threads - 1) / threads, ndm);
    size_t shmem = (size_t)nchan * sizeof(int);
    if (shmem > 48 * 1024) CK(cudaFuncSetAttribute(dedisperse_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)shmem));
    dedisperse_kernel<<<g, threads, shmem>>>(d_tr.p, nsamp, nchan, d_del.p, killmask ? d_kill.p : nullptr, nout, d_out.p);
    CK(cudaGetLastError());
    CK(cudaMemcpy(out, d_out.p, (size_t)ndm * nout * sizeof(float), cudaMemcpyDeviceToHost));
    return 0;
}

extern "C" int fbcuda_zero_dm(float *data, int nsamp, int nchan, int keep_channel_mean) {
    DevBuf<float> d, mean;
    size_t n = (size_t)nsamp * nchan;
    CK(d.alloc(n));
    CK(cudaMemcpy(d.p, data, n * sizeof(float), cudaMemcpyHostToDevice));
    if (keep_channel_mean) {
        CK(mean.alloc(nchan));
        channel_mean_kernel<<<(nchan + 127) / 128, 128>>>(d.p, nsamp, nchan, mean.p);
        CK(cudaGetLastError());
    }
    zero_dm_kernel<<<nsamp, 256>>>(d.p, nsamp, nchan, keep_channel_mean ? mean.p : nullptr);
    CK(cudaGetLastError());
    CK(cudaMemcpy(data, d.p, n * sizeof(float), cudaMemcpyDeviceToHost));
    return 0;
}

static void boxcar_noise(const double *cs, int n, int w, double *off, double *norm) {
    int m = n + 1 - w;
    int step = m / 16384 > 1 ? m / 16384 : 1;
    std::vector<double> y;
    for (int i = 0; i < m; i += step) y.push_back(cs[i + w] - cs[i]);
    *off = 0.0;
    *norm = 1.0 / std::sqrt((double)w);
    if (y.size() < 32) return;
    size_t mid = y.size() / 2;
    std::nth_element(y.begin(), y.begin() + mid, y.end());
    double med = y[mid];
    for (auto &v : y) v = std::fabs(v - med);
    std::nth_element(y.begin(), y.begin() + mid, y.end());
    double sd = 1.4826 * y[mid];
    *off = med;
    if (sd > 0.0) *norm = 1.0 / sd;
}

extern "C" int fbcuda_boxcar_snr(const float *series, int ndm, int nout, const int *widths, int nwidth,
                                 float *best_snr, int *best_width) {
    size_t rows = (size_t)ndm;
    std::string tmp;
    double *cs_h = nullptr;
    CK(cudaMallocHost(&cs_h, rows * (nout + 1) * sizeof(double)));
    for (size_t d = 0; d < rows; ++d) {
        double acc = 0.0;
        cs_h[d * (nout + 1)] = 0.0;
        for (int t = 0; t < nout; ++t) {
            acc += series[d * nout + t];
            cs_h[d * (nout + 1) + t + 1] = acc;
        }
    }
    std::vector<double> noise_h(rows * nwidth * 2);
    for (size_t d = 0; d < rows; ++d)
        for (int k = 0; k < nwidth; ++k)
            if (widths[k] <= nout) boxcar_noise(cs_h + d * (nout + 1), nout, widths[k], &noise_h[(d * nwidth + k) * 2], &noise_h[(d * nwidth + k) * 2 + 1]);
    DevBuf<double> cs, nz;
    DevBuf<int> w;
    DevBuf<float> bs;
    DevBuf<int> bw;
    cudaError_t e = cs.alloc(rows * (nout + 1));
    if (e != cudaSuccess) {
        cudaFreeHost(cs_h);
        CK(e);
    }
    CK(cudaMemcpy(cs.p, cs_h, rows * (nout + 1) * sizeof(double), cudaMemcpyHostToDevice));
    cudaFreeHost(cs_h);
    CK(nz.alloc(noise_h.size()));
    CK(cudaMemcpy(nz.p, noise_h.data(), noise_h.size() * sizeof(double), cudaMemcpyHostToDevice));
    CK(w.alloc(nwidth));
    CK(bs.alloc(rows * nout));
    CK(bw.alloc(rows * nout));
    CK(cudaMemcpy(w.p, widths, nwidth * sizeof(int), cudaMemcpyHostToDevice));
    dim3 g((nout + 255) / 256, ndm);
    boxcar_kernel<<<g, 256>>>(cs.p, ndm, nout, w.p, nwidth, nz.p, bs.p, bw.p);
    CK(cudaGetLastError());
    CK(cudaMemcpy(best_snr, bs.p, rows * nout * sizeof(float), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(best_width, bw.p, rows * nout * sizeof(int), cudaMemcpyDeviceToHost));
    return 0;
}

__global__ void fdmt_init_kernel(const float *tr, int nsamp, const int *order, int nf, int ndt0, float *buf) {
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= nsamp) return;
    for (int row = blockIdx.y; row < nf; row += gridDim.y) {
        int c = order[row];
        float acc = 0.0f;
        for (int k = 0; k < ndt0; ++k) {
            if (c >= 0 && t >= k) acc += tr[(size_t)c * nsamp + t - k];
            buf[((size_t)row * ndt0 + k) * nsamp + t] = t >= k ? acc : 0.0f;
        }
    }
}

__global__ void fdmt_level_kernel(const float *src, float *dst, const int *lo, const int *hi, const int *shift, int nrows, int nt) {
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= nt) return;
    for (int r = blockIdx.y; r < nrows; r += gridDim.y) {
        float x = src[(size_t)lo[r] * nt + t];
        int s = shift[r];
        dst[(size_t)r * nt + t] = t >= s ? x + src[(size_t)hi[r] * nt + t - s] : x;
    }
}

__global__ void fdmt_gather_kernel(const float *dmt, int nt, const int *rows, const int *shifts, int ndm, int nout, float *out) {
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= nout) return;
    for (int d = blockIdx.y; d < ndm; d += gridDim.y) {
        int src = t + shifts[d];
        out[(size_t)d * nout + t] = src < nt ? dmt[(size_t)rows[d] * nt + src] : 0.0f;
    }
}

extern "C" int fbcuda_fdmt(const float *data, int nsamp, int nchan, const int *order, int nf, int ndt0, int nlevels,
                           const int *level_rows, const int *lo, const int *hi, const int *shift, int max_rows, int ndm,
                           const int *out_rows, const int *out_shift, int nout, float *out) {
    size_t total = 0;
    for (int l = 0; l < nlevels; ++l) total += level_rows[l];
    DevBuf<float> d_in, d_tr, a, b, d_out;
    DevBuf<int> d_order, d_lo, d_hi, d_shift, d_rows, d_oshift;
    size_t n = (size_t)nsamp * nchan;
    CK(d_in.alloc(n));
    CK(d_tr.alloc(n));
    CK(a.alloc((size_t)max_rows * nsamp));
    CK(b.alloc((size_t)max_rows * nsamp));
    CK(d_out.alloc((size_t)ndm * nout));
    CK(d_order.alloc(nf));
    CK(d_lo.alloc(total ? total : 1));
    CK(d_hi.alloc(total ? total : 1));
    CK(d_shift.alloc(total ? total : 1));
    CK(d_rows.alloc(ndm));
    CK(d_oshift.alloc(ndm));
    CK(cudaMemcpy(d_in.p, data, n * sizeof(float), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_order.p, order, nf * sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_lo.p, lo, total * sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_hi.p, hi, total * sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_shift.p, shift, total * sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_rows.p, out_rows, ndm * sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_oshift.p, out_shift, ndm * sizeof(int), cudaMemcpyHostToDevice));
    dim3 tb(32, 8), tg((nchan + 31) / 32, (nsamp + 31) / 32);
    transpose_kernel<<<tg, tb>>>(d_in.p, d_tr.p, nsamp, nchan);
    CK(cudaGetLastError());
    int threads = 256;
    int gx = (nsamp + threads - 1) / threads;
    fdmt_init_kernel<<<dim3(gx, nf < 65535 ? nf : 65535), threads>>>(d_tr.p, nsamp, d_order.p, nf, ndt0, a.p);
    CK(cudaGetLastError());
    float *src = a.p, *dst = b.p;
    size_t off = 0;
    for (int l = 0; l < nlevels; ++l) {
        int nr = level_rows[l];
        fdmt_level_kernel<<<dim3(gx, nr < 65535 ? nr : 65535), threads>>>(src, dst, d_lo.p + off, d_hi.p + off, d_shift.p + off, nr, nsamp);
        CK(cudaGetLastError());
        off += nr;
        float *t = src;
        src = dst;
        dst = t;
    }
    fdmt_gather_kernel<<<dim3((nout + threads - 1) / threads, ndm < 65535 ? ndm : 65535), threads>>>(src, nsamp, d_rows.p, d_oshift.p, ndm, nout, d_out.p);
    CK(cudaGetLastError());
    CK(cudaMemcpy(out, d_out.p, (size_t)ndm * nout * sizeof(float), cudaMemcpyDeviceToHost));
    return 0;
}

#define CF(x)                                                       \
    do {                                                            \
        if ((x) != CUFFT_SUCCESS) {                                 \
            g_err = std::string(#x) + ": cuFFT failure";            \
            return -1;                                              \
        }                                                           \
    } while (0)

__global__ void row_demean_kernel(float *x, int n) {
    __shared__ double part[256];
    float *row = x + (size_t)blockIdx.x * n;
    double s = 0.0;
    for (int i = threadIdx.x; i < n; i += blockDim.x) s += row[i];
    part[threadIdx.x] = s;
    __syncthreads();
    for (int k = blockDim.x / 2; k > 0; k >>= 1) {
        if (threadIdx.x < k) part[threadIdx.x] += part[threadIdx.x + k];
        __syncthreads();
    }
    float m = (float)(part[0] / n);
    for (int i = threadIdx.x; i < n; i += blockDim.x) row[i] -= m;
}

__global__ void dered_block_kernel(cufftComplex *spec, int nbins) {
    __shared__ float v[256];
    cufftComplex *row = spec + (size_t)blockIdx.y * nbins;
    int lo = blockIdx.x * 256;
    int i = lo + threadIdx.x;
    v[threadIdx.x] = (i < nbins && i > 0) ? row[i].x * row[i].x + row[i].y * row[i].y : 3.0e38f;
    __syncthreads();
    for (int k = 2; k <= 256; k <<= 1) {
        for (int j = k >> 1; j > 0; j >>= 1) {
            int p = threadIdx.x ^ j;
            if (p > (int)threadIdx.x) {
                bool up = (threadIdx.x & k) == 0;
                float a = v[threadIdx.x], b = v[p];
                if ((a > b) == up) {
                    v[threadIdx.x] = b;
                    v[p] = a;
                }
            }
            __syncthreads();
        }
    }
    int valid = nbins - lo;
    if (valid > 256) valid = 256;
    if (lo == 0) valid -= 1;
    float med = valid > 0 ? v[valid / 2] : 1.0f;
    float scale = rsqrtf(fmaxf(med / 0.69314718f, 1e-30f));
    if (i < nbins) {
        if (i == 0) {
            row[i].x = 0.0f;
            row[i].y = 0.0f;
        } else {
            row[i].x *= scale;
            row[i].y *= scale;
        }
    }
}

__global__ void interbin_kernel(const cufftComplex *spec, int nbins, float *ib) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nbins) return;
    const cufftComplex *row = spec + (size_t)blockIdx.y * nbins;
    float *o = ib + (size_t)blockIdx.y * 2 * nbins;
    cufftComplex a = row[i];
    o[2 * i] = a.x * a.x + a.y * a.y;
    if (i + 1 < nbins) {
        float dx = 0.78539816f * (a.x - row[i + 1].x), dy = 0.78539816f * (a.y - row[i + 1].y);
        o[2 * i + 1] = dx * dx + dy * dy;
    } else {
        o[2 * i + 1] = 0.0f;
    }
}

__global__ void harmsum_kernel(const float *ib, int nib, int h, int stage, int lo, int hi, float thr, int cap,
                               int *count, int *c_row, int *c_stage, int *c_bin, float *c_pow) {
    int i = lo + blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= hi) return;
    const float *row = ib + (size_t)blockIdx.y * nib;
    float acc = 0.0f;
    for (int k = 1; k <= h; ++k) {
        long src = lroundf((float)((double)i * k / h));
        if (src < nib) acc += row[src];
    }
    if (acc >= thr) {
        int slot = atomicAdd(count, 1);
        if (slot < cap) {
            c_row[slot] = blockIdx.y;
            c_stage[slot] = stage;
            c_bin[slot] = i;
            c_pow[slot] = acc;
        }
    }
}

extern "C" int fbcuda_fft_search(const float *series, int nrows, int n, int nstages, const int *harms, const int *lo,
                                 const int *hi, const float *thr, int cap, int *c_row, int *c_stage, int *c_bin,
                                 float *c_pow) {
    int nbins = n / 2 + 1;
    DevBuf<float> d_x, d_ib, d_pow;
    DevBuf<cufftComplex> d_spec;
    DevBuf<int> d_count, d_row, d_stage, d_bin;
    CK(d_x.alloc((size_t)nrows * n));
    CK(d_spec.alloc((size_t)nrows * nbins));
    CK(d_ib.alloc((size_t)nrows * 2 * nbins));
    CK(d_count.alloc(1));
    CK(d_row.alloc(cap));
    CK(d_stage.alloc(cap));
    CK(d_bin.alloc(cap));
    CK(d_pow.alloc(cap));
    CK(cudaMemset(d_count.p, 0, sizeof(int)));
    CK(cudaMemcpy(d_x.p, series, (size_t)nrows * n * sizeof(float), cudaMemcpyHostToDevice));
    row_demean_kernel<<<nrows, 256>>>(d_x.p, n);
    CK(cudaGetLastError());
    cufftHandle plan;
    CF(cufftPlan1d(&plan, n, CUFFT_R2C, nrows));
    cufftResult fr = cufftExecR2C(plan, d_x.p, d_spec.p);
    cufftDestroy(plan);
    CF(fr);
    dered_block_kernel<<<dim3((nbins + 255) / 256, nrows), 256>>>(d_spec.p, nbins);
    CK(cudaGetLastError());
    interbin_kernel<<<dim3((nbins + 255) / 256, nrows), 256>>>(d_spec.p, nbins, d_ib.p);
    CK(cudaGetLastError());
    for (int s = 0; s < nstages; ++s) {
        int span = hi[s] - lo[s];
        if (span <= 0) continue;
        harmsum_kernel<<<dim3((span + 255) / 256, nrows), 256>>>(d_ib.p, 2 * nbins, harms[s], s, lo[s], hi[s], thr[s], cap,
                                                                 d_count.p, d_row.p, d_stage.p, d_bin.p, d_pow.p);
        CK(cudaGetLastError());
    }
    int count = 0;
    CK(cudaMemcpy(&count, d_count.p, sizeof(int), cudaMemcpyDeviceToHost));
    int m = count < cap ? count : cap;
    CK(cudaMemcpy(c_row, d_row.p, m * sizeof(int), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(c_stage, d_stage.p, m * sizeof(int), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(c_bin, d_bin.p, m * sizeof(int), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(c_pow, d_pow.p, m * sizeof(float), cudaMemcpyDeviceToHost));
    return count;
}
