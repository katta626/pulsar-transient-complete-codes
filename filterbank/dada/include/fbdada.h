#ifndef FBDADA_H
#define FBDADA_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int fbdada_create(int key, uint64_t nbufs, uint64_t bufsz, uint64_t hdr_nbufs, uint64_t hdr_bufsz, unsigned nreaders);
int fbdada_destroy(int key);
void *fbdada_open_reader(int key);
void *fbdada_open_writer(int key);
const char *fbdada_next_header(void *handle);
int64_t fbdada_read(void *handle, char *buf, uint64_t nbytes);
int fbdada_write_header(void *handle, const char *text);
int64_t fbdada_write(void *handle, const char *buf, uint64_t nbytes);
uint64_t fbdada_bufsz(void *handle);
void fbdada_close(void *handle);

#ifdef __cplusplus
}
#endif

#endif
