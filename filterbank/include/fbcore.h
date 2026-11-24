#ifndef FBCORE_H
#define FBCORE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Handle fb_handle;

const char *fb_version(void);
const char *fb_last_error(void);
void fb_free_string(char *p);

fb_handle *fb_open(const char *paths_json, const char *opts_json);
void fb_close(fb_handle *h);
char *fb_header(fb_handle *h);
int64_t fb_read(fb_handle *h, uint64_t start, uint64_t nsamp, float *out, uint64_t cap);

int64_t fb_write_fil(const char *path, const char *header_json, uint32_t nbits, const float *data, uint64_t nvalues);

char *fb_plan(const char *header_json, double lo, double hi, double step, double width, double tol);
char *fb_plan_from_list(const char *header_json, const char *dms_json);
int64_t fb_delay_table(const char *plan_json, int32_t *out, uint64_t cap);
int64_t fb_dedisperse(const char *header_json, const char *plan_json, const float *data, uint64_t nsamp,
                      const uint8_t *killmask, uint32_t nsub, float *out, uint64_t cap);

int64_t fb_iqrm(const double *x, uint64_t n, uint32_t radius, double threshold, uint8_t *out);
char *fb_mask_new(const char *header_json, double interval);
char *fb_rfi_block(const char *header_json, float *data, uint64_t nsamp, uint64_t start_sample,
                   const char *mask_json, const char *cfg_json);
int64_t fb_mask_apply(const char *header_json, float *data, uint64_t nsamp, uint64_t start_sample,
                      const char *mask_json, uint8_t fill_mean);
char *fb_presto_mask_read(const char *path, const char *header_json);
int64_t fb_presto_mask_write(const char *mask_json, double tsamp, const char *path);

char *fb_dat_info(const char *path);
char *fb_fdmt_setup(const char *header_json, const char *plan_json, uint64_t nsamp, const uint8_t *killmask);
double fb_sigma_from_power(double power, uint32_t numharm, double ntrials);
char *fb_dada_layout(const char *header_text, int sum_pols);
int64_t fb_dada_decode(const char *header_text, int sum_pols, const uint8_t *bytes, uint64_t nbytes, float *out, uint64_t cap);
int64_t fb_dat_read(const char *path, float *out, uint64_t cap);
int64_t fb_dat_write(const char *root, const char *header_json, double dm, double f_ref, const float *data, uint64_t n);

char *fb_sp_search(const char *header_json, const char *plan_json, const float *trials, uint64_t ndm, uint64_t nout,
                   uint64_t start_sample, const char *params_json);
char *fb_run_sp(const char *paths_json, const char *cfg_json);

#ifdef __cplusplus
}
#endif

#endif
