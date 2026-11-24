#ifndef FBCUDA_H
#define FBCUDA_H

#ifdef __cplusplus
extern "C" {
#endif

int fbcuda_device_count(void);
const char *fbcuda_last_error(void);
int fbcuda_set_device(int device);
int fbcuda_dedisperse(const float *data, int nsamp, int nchan, const int *delays, int ndm,
                      const unsigned char *killmask, int nout, float *out);
int fbcuda_zero_dm(float *data, int nsamp, int nchan, int keep_channel_mean);
int fbcuda_boxcar_snr(const float *series, int ndm, int nout, const int *widths, int nwidth,
                      float *best_snr, int *best_width);
int fbcuda_fdmt(const float *data, int nsamp, int nchan, const int *order, int nf, int ndt0, int nlevels,
                const int *level_rows, const int *lo, const int *hi, const int *shift, int max_rows, int ndm,
                const int *out_rows, const int *out_shift, int nout, float *out);
int fbcuda_fft_search(const float *series, int nrows, int n, int nstages, const int *harms, const int *lo,
                      const int *hi, const float *thr, int cap, int *c_row, int *c_stage, int *c_bin, float *c_pow);

#ifdef __cplusplus
}
#endif

#endif
