#include <cuda_runtime.h>

#include <cstdio>
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

__global__ void boxcar_kernel(const double *cs, int ndm, int nout, const int *widths, int nwidth,
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
        float s = (float)((row[t + w] - row[t]) * rsqrt((double)w));
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
    DevBuf<double> cs;
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
    CK(w.alloc(nwidth));
    CK(bs.alloc(rows * nout));
    CK(bw.alloc(rows * nout));
    CK(cudaMemcpy(w.p, widths, nwidth * sizeof(int), cudaMemcpyHostToDevice));
    dim3 g((nout + 255) / 256, ndm);
    boxcar_kernel<<<g, 256>>>(cs.p, ndm, nout, w.p, nwidth, bs.p, bw.p);
    CK(cudaGetLastError());
    CK(cudaMemcpy(best_snr, bs.p, rows * nout * sizeof(float), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(best_width, bw.p, rows * nout * sizeof(int), cudaMemcpyDeviceToHost));
    return 0;
}
