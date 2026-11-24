#pragma once
#include <algorithm>
#include <atomic>
#include <barrier>
#include <cmath>
#include <complex>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

struct dim3 {
    unsigned x, y, z;
    dim3(unsigned a = 1, unsigned b = 1, unsigned c = 1) : x(a), y(b), z(c) {}
};
struct uint3 {
    unsigned x, y, z;
};

inline thread_local uint3 threadIdx;
inline thread_local uint3 blockIdx;
inline thread_local dim3 blockDim;
inline thread_local dim3 gridDim;

#define __global__
#define __device__
#define __host__
#define __forceinline__ inline
#define __restrict__
#define __shared__ static

typedef int cudaError_t;
enum { cudaSuccess = 0, cudaErrorMemoryAllocation = 2 };
enum cudaMemcpyKind { cudaMemcpyHostToHost, cudaMemcpyHostToDevice, cudaMemcpyDeviceToHost, cudaMemcpyDeviceToDevice };
enum cudaFuncAttribute { cudaFuncAttributeMaxDynamicSharedMemorySize };

template <class T>
inline cudaError_t cudaMalloc(T **p, size_t n) {
    *p = (T *)std::calloc(n ? n : 1, 1);
    return *p ? cudaSuccess : cudaErrorMemoryAllocation;
}
template <class T>
inline cudaError_t cudaMallocHost(T **p, size_t n) {
    return cudaMalloc(p, n);
}
inline cudaError_t cudaFree(void *p) {
    std::free(p);
    return cudaSuccess;
}
inline cudaError_t cudaFreeHost(void *p) { return cudaFree(p); }
inline cudaError_t cudaMemcpy(void *d, const void *s, size_t n, cudaMemcpyKind) {
    std::memcpy(d, s, n);
    return cudaSuccess;
}
inline cudaError_t cudaMemset(void *d, int v, size_t n) {
    std::memset(d, v, n);
    return cudaSuccess;
}
inline cudaError_t cudaGetLastError() { return cudaSuccess; }
inline cudaError_t cudaDeviceSynchronize() { return cudaSuccess; }
inline const char *cudaGetErrorString(cudaError_t) { return "emulated cuda error"; }
inline cudaError_t cudaGetDeviceCount(int *n) {
    *n = 1;
    return cudaSuccess;
}
inline cudaError_t cudaSetDevice(int) { return cudaSuccess; }
template <class F>
inline cudaError_t cudaFuncSetAttribute(F, cudaFuncAttribute, int) {
    return cudaSuccess;
}

inline float rsqrtf(float x) { return 1.0f / std::sqrt(x); }
inline double rsqrt(double x) { return 1.0 / std::sqrt(x); }
inline float atomicAdd(float *p, float v) { return std::atomic_ref<float>(*p).fetch_add(v); }
inline int atomicAdd(int *p, int v) { return std::atomic_ref<int>(*p).fetch_add(v); }
inline unsigned atomicAdd(unsigned *p, unsigned v) { return std::atomic_ref<unsigned>(*p).fetch_add(v); }

struct FbemuBlockSync {
    std::mutex m;
    std::condition_variable cv;
    int active = 0, waiting = 0;
    long gen = 0;
    void sync() {
        std::unique_lock<std::mutex> l(m);
        long g = gen;
        if (++waiting >= active) {
            waiting = 0;
            ++gen;
            cv.notify_all();
        } else {
            cv.wait(l, [&] { return gen != g; });
        }
    }
    void exit() {
        std::unique_lock<std::mutex> l(m);
        --active;
        if (active > 0 && waiting >= active) {
            waiting = 0;
            ++gen;
            cv.notify_all();
        }
    }
};

inline FbemuBlockSync *fbemu_sync = nullptr;
inline std::vector<unsigned char> fbemu_dyn;
inline void __syncthreads() { fbemu_sync->sync(); }
inline void *fbemu_dynamic_shared() { return fbemu_dyn.data(); }

struct FbemuCfg {
    dim3 g, b;
    size_t sh;
};
inline FbemuCfg fbemu_cfg(dim3 g, dim3 b, size_t sh = 0, int = 0) { return {g, b, sh}; }

template <class K, class... A>
inline void fbemu_launch(FbemuCfg c, K k, A... args) {
    static std::mutex launch_lock;
    std::lock_guard<std::mutex> guard(launch_lock);
    size_t nt = (size_t)c.b.x * c.b.y * c.b.z;
    size_t nb = (size_t)c.g.x * c.g.y * c.g.z;
    fbemu_dyn.assign(c.sh ? c.sh : 16, 0);
    FbemuBlockSync bs;
    bs.active = (int)nt;
    fbemu_sync = &bs;
    auto reset = [&]() noexcept {
        bs.active = (int)nt;
        bs.waiting = 0;
    };
    std::barrier end((std::ptrdiff_t)nt, reset);
    std::vector<std::thread> pool;
    for (size_t i = 0; i < nt; ++i) {
        pool.emplace_back([&, i] {
            blockDim = c.b;
            gridDim = c.g;
            threadIdx = {(unsigned)(i % c.b.x), (unsigned)((i / c.b.x) % c.b.y), (unsigned)(i / ((size_t)c.b.x * c.b.y))};
            for (size_t blk = 0; blk < nb; ++blk) {
                blockIdx = {(unsigned)(blk % c.g.x), (unsigned)((blk / c.g.x) % c.g.y), (unsigned)(blk / ((size_t)c.g.x * c.g.y))};
                k(args...);
                bs.exit();
                end.arrive_and_wait();
            }
        });
    }
    for (auto &t : pool) t.join();
    fbemu_sync = nullptr;
}
