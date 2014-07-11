/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013, 2014 Sebastian Wouters

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

#ifndef CHEMPS2_OPTIONS
#define CHEMPS2_OPTIONS

#include <stdlib.h>
#include <string>

using std::string;

namespace CheMPS2{

   const bool   DMRGSCF_debugPrint            = false;
   const bool   DMRGSCF_storeUnitary          = true;
   const bool   DMRGSCF_rotate2DMtoNO         = true;
   const double DMRGSCF_gradientNormThreshold = 1e-8;
   const string DMRGSCF_unitaryStorageName    = "CheMPS2_CASSCF.h5";
   const int    DMRGSCF_maxlinsizeCutoff      = 100;
   
   const string DMRGSCF_DIISstorageName       = "CheMPS2_DIIS.h5";
   const bool   DMRGSCF_storeDIIS             = true;
   const double DMRGSCF_DIISgradientBranch    = 5e-3;
   const int    DMRGSCF_numDIISvecs           = 10;

   const string TMPpath                       = "/tmp/";
   const bool   DMRG_printDiscardedWeight     = false;
   const bool   DMRG_storeRenormOptrOnDisk    = true;
   const bool   DMRG_storeMpsOnDisk           = false;
   
   const bool   HAMILTONIAN_debugPrint        = false;
   const string HAMILTONIAN_TmatStorageName   = "CheMPS2_Ham_Tmat.h5";
   const string HAMILTONIAN_VmatStorageName   = "CheMPS2_Ham_Vmat.h5"; 
   const string HAMILTONIAN_ParentStorageName = "CheMPS2_Ham_parent.h5";
   
   const string TWODM_2DM_A_storagename       = "CheMPS2_2DM-A.h5";
   const string TWODM_2DM_B_storagename       = "CheMPS2_2DM-B.h5";
   
   const bool   HEFF_debugPrint               = true;
   const int    HEFF_DAVIDSON_NUM_VEC         = 32;
   const int    HEFF_DAVIDSON_NUM_VEC_KEEP    = 3;
   const double HEFF_DAVIDSON_PRECOND_CUTOFF  = 1e-12;
   const double HEFF_DAVIDSON_RTOL_BASE       = 1e-10;
   
   const bool   SYBK_debugPrint               = false;
   const int    SYBK_dimensionCutoff          = 262144;
   
   const double TENSORT_orthoComparison       = 1e-13;

}

#endif

