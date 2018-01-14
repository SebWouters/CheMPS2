/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2018 Sebastian Wouters

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

#ifndef OPTIONS_CHEMPS2_H
#define OPTIONS_CHEMPS2_H

#include <stdlib.h>
#include <string>

using std::string;

namespace CheMPS2{

   const int    DMRGSCF_maxIterations         = 100;
   const double DMRGSCF_gradientNormThreshold = 1e-6;
   const bool   DMRGSCF_storeUnitary          = true;
   const string DMRGSCF_unitary_storage_name  = "CheMPS2_CASSCF.h5";
   const string DMRGSCF_eri_storage_name      = "CheMPS2_eri_temp.h5";
   const string DMRGSCF_f4rdm_name            = "CheMPS2_f4rdm.h5";
   const int    DMRGSCF_max_mem_eri_tfo       = 100 * 100 * 100 * 100; // Measured in number of doubles
   const bool   DMRGSCF_debugPrint            = false;
   const bool   DMRGSCF_stateAveraged         = true;

   const int    DMRGSCF_whichActiveSpace      = 0;
   const bool   DMRGSCF_dumpCorrelations      = false;
   const bool   DMRGSCF_startLocRandom        = false;

   const bool   DMRGSCF_doDIIS                = false;
   const double DMRGSCF_DIISgradientBranch    = 1e-2;
   const int    DMRGSCF_numDIISvecs           = 7;
   const bool   DMRGSCF_storeDIIS             = true;
   const string DMRGSCF_diis_storage_name     = "CheMPS2_DIIS.h5";

   const double CASPT2_OVLP_CUTOFF            = 1e-8;

   const double CONJ_GRADIENT_RTOL            = 1e-10;
   const double CONJ_GRADIENT_PRECOND_CUTOFF  = 1e-12;

   const string defaultTMPpath                = "/tmp";
   const bool   DMRG_storeRenormOptrOnDisk    = true;
   const bool   DMRG_storeMpsOnDisk           = false;
   const string DMRG_MPS_storage_prefix       = "CheMPS2_MPS";
   const string DMRG_OPERATOR_storage_prefix  = "CheMPS2_Operators_";

   const bool   HAMILTONIAN_debugPrint        = false;
   const string HAMILTONIAN_TmatStorageName   = "CheMPS2_Ham_Tmat.h5";
   const string HAMILTONIAN_VmatStorageName   = "CheMPS2_Ham_Vmat.h5";
   const string HAMILTONIAN_ParentStorageName = "CheMPS2_Ham_parent.h5";

   const string TWO_RDM_storagename           = "CheMPS2_2DM.h5";
   const string THREE_RDM_storage_prefix      = "CheMPS2_3DM_";

   const bool   HEFF_debugPrint               = true;
   const int    DAVIDSON_NUM_VEC              = 32;
   const int    DAVIDSON_NUM_VEC_KEEP         = 3;
   const double DAVIDSON_PRECOND_CUTOFF       = 1e-12;
   const double DAVIDSON_FCI_RTOL             = 1e-10;  // Base value for FCI and augmented Hessian diagonalization
   const double DAVIDSON_DMRG_RTOL            = 1e-5;   // Block's Davidson tolerance would correspond to HEFF_DAVIDSON_DMRG_RTOL^2

   const int    SYBK_dimensionCutoff          = 262144;

   const double TENSORT_orthoComparison       = 1e-13;

   const bool   CORRELATIONS_debugPrint       = false;
   const double CORRELATIONS_discardEig       = 1e-100;

   const double EDMISTONRUED_gradThreshold    = 1e-8;
   const int    EDMISTONRUED_maxIter          = 1000;
   const int    EDMISTONRUED_maxIterBackTfo   = 15;

}

#endif

