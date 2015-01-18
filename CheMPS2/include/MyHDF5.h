/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2015 Sebastian Wouters

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef HDF5_CHEMPS2_H
#define HDF5_CHEMPS2_H

   //Force the use of the 1.8 API of HDF5
   #undef H5_USE_16_API
   #define H5_NO_DEPRECATED_SYMBOLS
   #define H5Acreate_vers 2
   #define H5Dcreate_vers 2
   #define H5Dopen_vers 2
   #define H5Gcreate_vers 2
   #define H5Gopen_vers 2

   #include <hdf5.h>

#endif /* HDF5_CHEMPS2_H */
