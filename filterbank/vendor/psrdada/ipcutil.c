#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "ipcutil.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>
#include <sys/mman.h>

#ifdef __APPLE__
static int semtimedop(int semid, struct sembuf *ops, size_t nops, const struct timespec *ts)
{
  struct sembuf b = *ops;
  b.sem_flg |= IPC_NOWAIT;
  uint64_t waited = 0;
  uint64_t limit = (uint64_t) ts->tv_sec * 1000000000ULL + (uint64_t) ts->tv_nsec;
  struct timespec step = {0, 100000};
  for (;;)
  {
    if (semop(semid, &b, nops) == 0)
      return 0;
    if (errno != EAGAIN || waited >= limit)
      return -1;
    nanosleep(&step, NULL);
    waited += 100000;
  }
}
#endif

// #define _DEBUG 1

/* *************************************************************** */
/*!
  Returns a shared memory block and its shmid
*/
void* ipc_alloc (key_t key, size_t size, int flag, int* shmid)
{
  void* buf = 0;
  int id = 0;
#ifdef _DEBUG
  fprintf (stderr, "ipc_alloc: key=%x size=%ld flag=%x\n", key, size, flag);
#endif

  id = shmget (key, size, flag);
  if (id < 0) {
    fprintf (stderr, "ipc_alloc: shmget (key=%x, size=%ld, flag=%x) %s\n",
             key, size, flag, strerror(errno));
    return 0;
  }

#ifdef _DEBUG
  fprintf (stderr, "ipc_alloc: shmid=%d\n", id);
#endif

  buf = shmat (id, 0, flag);

  if (buf == (void*)-1) {
    fprintf (stderr,
	     "ipc_alloc: shmat (shmid=%d) %s\n"
	     "ipc_alloc: after shmget (key=%x, size=%ld, flag=%x)\n",
	     id, strerror(errno), key, size, flag);
    return 0;
  }

#ifdef _DEBUG
  fprintf (stderr, "ipc_alloc: shmat=%p\n", buf);
#endif

  if (shmid)
    *shmid = id;

  return buf;
}

int ipc_semop (int semid, short num, short op, short flag)
{
  struct sembuf semopbuf;

  semopbuf.sem_num = num;
  semopbuf.sem_op = op;
  semopbuf.sem_flg = flag;

  if (semop (semid, &semopbuf, 1) < 0) {
    if (!(flag | IPC_NOWAIT))
      perror ("ipc_semop: semop");
    return -1;
  }
  return 0;
}

int ipc_semtimedop (int semid, short num, short op, short flag, uint64_t timeout)
{
  if (timeout == 0)
  {
    return ipc_semop (semid, num, op, flag);
  }

  struct sembuf semopbuf;
  struct timespec ts;
  const uint64_t ns_per_second = 1000000000;

  ts.tv_sec = 0;
  if (timeout >= ns_per_second)
  {
    ts.tv_sec = timeout / ns_per_second;
    ts.tv_nsec = timeout - (ts.tv_sec * ns_per_second);
  }
  else
  {
    ts.tv_nsec = timeout;
  }

  semopbuf.sem_num = num;
  semopbuf.sem_op = op;
  semopbuf.sem_flg = flag;

  if (semtimedop (semid, &semopbuf, 1, &ts) < 0)
  {
    if ( errno == EAGAIN )
    {
      return 1;
    }
    if (!(flag | IPC_NOWAIT))
    {
      perror ("ipc_semop: semop");
    }
    return -1;
  }
  return 0;
}

int validate_shm_create_flags(int shm_create_flags)
{
  // Allocate a one-page shared memory region with the provided flags to ensure nothing errors out
  const size_t page_size = getpagesize();
  int id = shmget(IPC_PRIVATE, page_size, shm_create_flags);
  if (id < 0)
  {
    perror("validate_shm_create_flags: shmget");
    return -1;
  }

  // Then clean up after ourselves
  if (shmctl(id, IPC_RMID, NULL) < 0)
  {
    perror("validate_shm_create_flags: shmctl(IPC_RMID)");
    return -1;
  }

  return 0;
}
