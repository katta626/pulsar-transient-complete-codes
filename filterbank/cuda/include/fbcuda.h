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

#ifdef __cplusplus
}
#endif

#endif
