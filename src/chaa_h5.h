/* chaa_h5.h — minimal prototypes for the HDF5 entry points used by
 * Chaa's writer (Hdf5IO.chpl).  Matches the HDF5 >= 1.10 ABI, where
 * hid_t is a 64-bit integer; avoids a build-time dependency on the
 * full hdf5.h header. */
#ifndef CHAA_H5_H
#define CHAA_H5_H

#include <stdint.h>

/* GPU builds compile the generated glue as C++ — keep C linkage */
#ifdef __cplusplus
extern "C" {
#endif

int     H5open(void);
int64_t H5Fcreate(const char *name, unsigned flags, int64_t fcpl,
                  int64_t fapl);
int     H5Fclose(int64_t f);
int64_t H5Screate_simple(int rank, const uint64_t *dims,
                         const uint64_t *maxdims);
int     H5Sclose(int64_t s);
int64_t H5Dcreate2(int64_t loc, const char *name, int64_t dtype,
                   int64_t space, int64_t lcpl, int64_t dcpl, int64_t dapl);
int     H5Dwrite(int64_t dset, int64_t memtype, int64_t memspace,
                 int64_t filespace, int64_t xfer, const void *buf);
int     H5Dclose(int64_t d);

extern int64_t H5T_NATIVE_DOUBLE_g;

#ifdef __cplusplus
}
#endif

#endif
