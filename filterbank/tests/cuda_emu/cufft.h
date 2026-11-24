#pragma once
#include <complex>
#include <vector>

#include "cuda_runtime.h"

typedef int cufftHandle;
typedef int cufftResult;
typedef float cufftReal;
struct cufftComplex {
    float x, y;
};
enum { CUFFT_SUCCESS = 0, CUFFT_R2C = 0x2a, CUFFT_C2C = 0x29, CUFFT_FORWARD = -1, CUFFT_INVERSE = 1 };

struct FbemuPlan {
    int n, type, batch;
};
inline std::vector<FbemuPlan> &fbemu_plans() {
    static std::vector<FbemuPlan> p;
    return p;
}
inline cufftResult cufftPlan1d(cufftHandle *h, int n, int type, int batch) {
    fbemu_plans().push_back({n, type, batch});
    *h = (int)fbemu_plans().size() - 1;
    return CUFFT_SUCCESS;
}
inline cufftResult cufftDestroy(cufftHandle) { return CUFFT_SUCCESS; }

inline void fbemu_fft(std::vector<std::complex<double>> &a, int sign) {
    size_t n = a.size();
    if (n & (n - 1)) {
        std::vector<std::complex<double>> o(n);
        for (size_t k = 0; k < n; ++k)
            for (size_t j = 0; j < n; ++j) o[k] += a[j] * std::polar(1.0, sign * 2 * M_PI * (double)(j * k % n) / n);
        a = o;
        return;
    }
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (size_t len = 2; len <= n; len <<= 1) {
        std::complex<double> w = std::polar(1.0, sign * 2 * M_PI / len);
        for (size_t i = 0; i < n; i += len) {
            std::complex<double> u = 1;
            for (size_t j = 0; j < len / 2; ++j) {
                auto x = a[i + j], y = a[i + j + len / 2] * u;
                a[i + j] = x + y;
                a[i + j + len / 2] = x - y;
                u *= w;
            }
        }
    }
}

inline cufftResult cufftExecR2C(cufftHandle h, cufftReal *in, cufftComplex *out) {
    auto p = fbemu_plans()[h];
    for (int b = 0; b < p.batch; ++b) {
        std::vector<std::complex<double>> a(p.n);
        for (int i = 0; i < p.n; ++i) a[i] = in[(size_t)b * p.n + i];
        fbemu_fft(a, -1);
        for (int i = 0; i <= p.n / 2; ++i) out[(size_t)b * (p.n / 2 + 1) + i] = {(float)a[i].real(), (float)a[i].imag()};
    }
    return CUFFT_SUCCESS;
}

inline cufftResult cufftExecC2C(cufftHandle h, cufftComplex *in, cufftComplex *out, int dir) {
    auto p = fbemu_plans()[h];
    for (int b = 0; b < p.batch; ++b) {
        std::vector<std::complex<double>> a(p.n);
        for (int i = 0; i < p.n; ++i) a[i] = {in[(size_t)b * p.n + i].x, in[(size_t)b * p.n + i].y};
        fbemu_fft(a, dir);
        for (int i = 0; i < p.n; ++i) out[(size_t)b * p.n + i] = {(float)a[i].real(), (float)a[i].imag()};
    }
    return CUFFT_SUCCESS;
}
