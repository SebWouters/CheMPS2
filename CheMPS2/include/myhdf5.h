#ifndef MY_HDF5_H
#define MY_HDF5_H

// force the use of the 1.8 API of HDF5
#undef H5_USE_16_API
#define H5_NO_DEPRECATED_SYMBOLS
#define H5Acreate_vers 2
#define H5Dcreate_vers 2
#define H5Dopen_vers 2
#define H5Gcreate_vers 2
#define H5Gopen_vers 2

#include <hdf5.h>

#endif /* MY_HDF5_H */
