#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fbdada.h"
#include "dada_hdu.h"
#include "ipcbuf.h"
#include "ipcio.h"
#include "multilog.h"

typedef struct {
  dada_hdu_t *hdu;
  multilog_t *log;
  int writer;
  int locked;
} fbdada_t;

int fbdada_create(int key, uint64_t nbufs, uint64_t bufsz, uint64_t hdr_nbufs, uint64_t hdr_bufsz, unsigned nreaders)
{
  ipcbuf_t data = IPCBUF_INIT;
  ipcbuf_t hdr = IPCBUF_INIT;
  if (ipcbuf_create(&data, key, nbufs, bufsz, nreaders) < 0)
    return -1;
  if (ipcbuf_create(&hdr, key + 1, hdr_nbufs, hdr_bufsz, nreaders) < 0)
  {
    ipcbuf_destroy(&data);
    return -2;
  }
  ipcbuf_disconnect(&data);
  ipcbuf_disconnect(&hdr);
  return 0;
}

int fbdada_destroy(int key)
{
  int rc = 0;
  ipcbuf_t data = IPCBUF_INIT;
  ipcbuf_t hdr = IPCBUF_INIT;
  if (ipcbuf_connect(&data, key) == 0)
    rc |= ipcbuf_destroy(&data);
  else
    rc = -1;
  if (ipcbuf_connect(&hdr, key + 1) == 0)
    rc |= ipcbuf_destroy(&hdr);
  else
    rc = -1;
  return rc;
}

static fbdada_t *connect_hdu(int key)
{
  fbdada_t *f = calloc(1, sizeof(fbdada_t));
  if (!f)
    return NULL;
  f->log = multilog_open("fbdada", 0);
  multilog_add(f->log, stderr);
  f->hdu = dada_hdu_create(f->log);
  dada_hdu_set_key(f->hdu, key);
  if (dada_hdu_connect(f->hdu) < 0)
  {
    dada_hdu_destroy(f->hdu);
    multilog_close(f->log);
    free(f);
    return NULL;
  }
  return f;
}

void *fbdada_open_reader(int key)
{
  fbdada_t *f = connect_hdu(key);
  if (!f)
    return NULL;
  if (dada_hdu_lock_read(f->hdu) < 0)
  {
    fbdada_close(f);
    return NULL;
  }
  f->locked = 1;
  return f;
}

void *fbdada_open_writer(int key)
{
  fbdada_t *f = connect_hdu(key);
  if (!f)
    return NULL;
  f->writer = 1;
  if (dada_hdu_lock_write(f->hdu) < 0)
  {
    fbdada_close(f);
    return NULL;
  }
  f->locked = 1;
  return f;
}

const char *fbdada_next_header(void *handle)
{
  fbdada_t *f = handle;
  if (f->hdu->header)
  {
    free(f->hdu->header);
    f->hdu->header = NULL;
  }
  if (dada_hdu_open(f->hdu) < 0)
    return NULL;
  return f->hdu->header;
}

int64_t fbdada_read(void *handle, char *buf, uint64_t nbytes)
{
  fbdada_t *f = handle;
  return ipcio_read(f->hdu->data_block, buf, nbytes);
}

int fbdada_write_header(void *handle, const char *text)
{
  fbdada_t *f = handle;
  uint64_t sz = ipcbuf_get_bufsz(f->hdu->header_block);
  char *h = ipcbuf_get_next_write(f->hdu->header_block);
  if (!h)
    return -1;
  memset(h, 0, sz);
  strncpy(h, text, sz - 1);
  return ipcbuf_mark_filled(f->hdu->header_block, sz);
}

int64_t fbdada_write(void *handle, const char *buf, uint64_t nbytes)
{
  fbdada_t *f = handle;
  return ipcio_write(f->hdu->data_block, (char *)buf, nbytes);
}

uint64_t fbdada_bufsz(void *handle)
{
  fbdada_t *f = handle;
  return ipcbuf_get_bufsz((ipcbuf_t *)f->hdu->data_block);
}

void fbdada_close(void *handle)
{
  fbdada_t *f = handle;
  if (!f)
    return;
  if (f->locked)
  {
    if (f->writer)
      dada_hdu_unlock_write(f->hdu);
    else
      dada_hdu_unlock_read(f->hdu);
  }
  dada_hdu_disconnect(f->hdu);
  dada_hdu_destroy(f->hdu);
  multilog_close(f->log);
  free(f);
}
