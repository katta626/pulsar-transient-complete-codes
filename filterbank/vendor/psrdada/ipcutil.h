#ifndef __DADA_IPCUTIL_H
#define __DADA_IPCUTIL_H

#include <sys/types.h>
#include <inttypes.h>

/*************************************************************************

   utilities for creation of shared memory and operations on semaphores

  ************************************************************************ */

#define IPCUTIL_PERM 0666 /* default: read/write permissions for all */

#define IPCUTIL_DEFAULT_CREATE_FLAGS (IPCUTIL_PERM | IPC_CREAT | IPC_EXCL) /* default flags for creating semaphores and shared memory regions */

#ifdef __cplusplus
extern "C" {
#endif

  /* allocate size bytes in shared memory with the specified flags and key.
     returns the pointer to the base address and the shmid, id */
  void* ipc_alloc (key_t key, size_t size, int flag, int* id);

  /* operate on the specified semaphore */
  int ipc_semop (int semid, short num, short incr, short flag);

  /**
   * @brief operate on the specified semaphore, within the timeout.
   *
   * If the timeout is zero, then the ipc_semop function will be used.
   *
   * @param semid the semaphore identifier on which to operated
   * @param num the semaphore number
   * @param op the semaphore operation
   * @param flag the semaphore operation flags (e.g. IPC_NOWAIT, SEM_UNDO)
   * @param timeout the maximum timeout in nano seconds
   * @return -1 on operation error, 0 if the semaphore operation has completed within the time, +1 if the timeout occurred
   */
  int ipc_semtimedop (int semid, short num, short op, short flag, uint64_t timeout);

  /**
   * @brief Validate the given flags to be passed to shmget(2) when creating a new shared memory segment, by creating and destroying
   * a 1-byte segment with those flags.
   *
   * @param shm_create_flags the flags to use when creating the shared memory segment.
   * @return 0 if the SHM creation and destruction was successful, -1 otherwise. On error, \c errno will be set to indicate the error. 
   */
  int validate_shm_create_flags(int shm_create_flags);

#ifdef __cplusplus
}
#endif

#endif // __DADA_IPCUTIL_H
