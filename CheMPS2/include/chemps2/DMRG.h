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

#ifndef DMRG_CHEMPS2_H
#define DMRG_CHEMPS2_H

#include <string>

#include "Options.h"
#include "Problem.h"
#include "TensorT.h"
#include "TensorX.h"
#include "TensorL.h"
#include "SyBookkeeper.h"
#include "TensorF1.h"
#include "TensorF0.h"
#include "TensorS0.h"
#include "TensorS1.h"
#include "TensorA.h"
#include "TensorB.h"
#include "TensorC.h"
#include "TensorD.h"
#include "TensorQ.h"
#include "TensorO.h"
#include "TwoDM.h"
#include "Correlations.h"
#include "Heff.h"
#include "Sobject.h"
#include "ConvergenceScheme.h"
#include "MyHDF5.h"

//For the timings of the different parts of DMRG
#define CHEMPS2_TIME_S_JOIN      0
#define CHEMPS2_TIME_S_SOLVE     1
#define CHEMPS2_TIME_S_SPLIT     2
#define CHEMPS2_TIME_TENS_TOTAL  3
#define CHEMPS2_TIME_TENS_ALLOC  4
#define CHEMPS2_TIME_TENS_FREE   5
#define CHEMPS2_TIME_TENS_DISKIO 6
#define CHEMPS2_TIME_TENS_CALC   7
#define CHEMPS2_TIME_VECLENGTH   8

namespace CheMPS2{
/** DMRG class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date July 31, 2013
    
    The DMRG class solves the Problem with its given parameters. A fully SU(2) symmetric MPS wavefunction is variationally optimized in a two-site sweep algorithm. When the solution has been reached, the converged energy and spin contracted 2DMs can be accessed. For more information, please take a look at
    
    S. Wouters, W. Poelmans, P.W. Ayers and D. Van Neck, \n
    CheMPS2: a free open-source spin-adapted implementation of the density matrix renormalization group for ab initio quantum chemistry, \n
    Computer Physics Communications 185, 1501-1514 (2014) \n
    http://dx.doi.org/10.1016/j.cpc.2014.01.019 \n
    http://arxiv.org/abs/1312.2415
    
    S. Wouters and D. Van Neck, \n
    The density matrix renormalization group for ab initio quantum chemistry, \n
    European Physical Journal D 68, 272 (2014) \n
    http://dx.doi.org/10.1140/epjd/e2014-50500-1 \n
    http://arxiv.org/abs/1407.2040
    
    The user manual: http://sebwouters.github.io/CheMPS2/index.html
    */
   class DMRG{

      public:
      
         //! Constructor
         /** \param Probin The problem to be solved
             \param OptSchemeIn The optimization scheme for the DMRG sweeps
             \param makechkpt Whether or not to save MPS checkpoints in the working directory
             \param tmpfolder Temporary folder on a large partition to store the renormalized operators on disk (by default "/tmp") */
         DMRG(Problem * Probin, ConvergenceScheme * OptSchemeIn, const bool makechkpt=CheMPS2::DMRG_storeMpsOnDisk, const string tmpfolder=CheMPS2::defaultTMPpath);
         
         //! Destructor
         virtual ~DMRG();
         
         //! Solver
         /** \return The min. energy encountered so far during the sweeps. */
         double Solve();
         
         //! Calculate the 2DM. Note that the DMRG class cannot be used for further updates anymore !!!
         void calc2DMandCorrelations();
         
         //! Get the pointer to the 2DM
         /** \return The 2DM. Returns a NULL pointer if not yet calculated. */
         TwoDM * get2DM(){ return the2DM; }
         
         //! Get the pointer to the Correlations
         /** \return The Correlations. Returns a NULL pointer if not yet calculated. */
         Correlations * getCorrelations(){ return theCorr; }
         
         //! Get a specific FCI coefficient. The array coeff contains the occupation numbers of the L Hamiltonian orbitals. It is assumed that the unpaired electrons are all alpha electrons, and that this number equals twice the total targeted spin.
         /** \param coeff Array containing the occupation numbers of the L Hamiltonian orbitals (occupations can be 0, 1, or 2).
             \return The desired FCI coefficient */
         double getSpecificCoefficient(int * coeff) const;
         
         //! Get a specific FCI coefficient. The arrays alpha and beta contain the alpha and beta occupation numbers of the L Hamiltonian orbitals.
         /** \param alpha Array containing the alpha electron occupation numbers of the L Hamiltonian orbitals (occupations can be 0 or 1).
             \param beta  Array containing the beta  electron occupation numbers of the L Hamiltonian orbitals (occupations can be 0 or 1).
             \param mpi_chemps2_master_only When running with MPI, whether only the master process should calculate the FCI coefficient. If false, any process can calculate the FCI coefficient, and it won't be broadcasted (allows for parallel FCI coefficient calculation).
             \return The desired FCI coefficient */
         double getFCIcoefficient(int * alpha, int * beta, const bool mpi_chemps2_master_only=true) const;
         
         //! Call "rm " + CheMPS2::DMRG_MPS_storage_prefix + "*.h5"
         void deleteStoredMPS();
         
         //! Call "rm " + tempfolder + "/" + CheMPS2::DMRG_OPERATOR_storage_prefix + string(thePID) + "*.h5";
         void deleteStoredOperators();
         
         //! Activate the necessary storage and machinery to handle excitations
         /** \param maxExcIn The max. number of excitations desired */
         void activateExcitations(const int maxExcIn);
         
         //! Push back current calculation and set everything up to calculate a (new) excitation
         /** \param EshiftIn To the Hamiltonian, a level shift is introduced to exclude the previously calculated MPS: Hnew = Hold + EshiftIn * | prev> <prev| */
         void newExcitation(const double EshiftIn);
         
         //! Print the license
         void PrintLicense();
         
      private:
      
         //Setup the DMRG SyBK and MPS (in separate function to allow pushbacks and recreations for excited states)
         void setupBookkeeperAndMPS();
      
         //! DMRG MPS + virt. dim. storage filename
         string MPSstoragename;
         
         //The optimization scheme for the DMRG sweeps (externally allocated, filled and deleted)
         ConvergenceScheme * OptScheme;
         
         //Whether or not the MPS was loaded from disk to memory at the start
         bool loadedMPS;
         
         //Integer to distinguish storage between different calculations
         int thePID;
      
         //Pointer to the Problem --> constructed and destructed outside of this class
         Problem * Prob;
         
         //The number of orbitals: copied here so that the DMRG destructor doesn't depend on whether Prob still exists
         int L;
         
         //Minimum energy encountered during all performed micro-iterations (as opposed to the 2DM/edge energy)
         double TotalMinEnergy;
         
         //Minimum energy encountered during the micro-iterations of the last performed sweep (as opposed to the 2DM/edge energy)
         double LastMinEnergy;
         
         //Max. discarded weight of last sweep
         double MaxDiscWeightLastSweep;
         
         //Symmetry information object
         SyBookkeeper * denBK;
         
         //The MPS
         TensorT ** MPS;
         
         //The TwoDM
         TwoDM * the2DM;
         
         //Whether the2DM is allocated
         bool the2DMallocated;
         
         //The Correlations
         Correlations * theCorr;
         
         //Whether the Correlations is allocated
         bool theCorrAllocated;
         
         //Whether or not allocated
         int * isAllocated;
         
         //TensorL's
         TensorL *** Ltensors;
         
         //TensorX's
         TensorX ** Xtensors;
         
         //TensorF0's
         TensorF0 **** F0tensors;
         
         //TensorF1's
         TensorF1 **** F1tensors;
         
         //TensorS0's
         TensorS0 **** S0tensors;
         
         //TensorS1's
         TensorS1 **** S1tensors;
         
         //TensorA's
         TensorA **** Atensors;
         
         //TensorB's
         TensorB **** Btensors;
         
         //TensorC's
         TensorC **** Ctensors;
         
         //TensorD's
         TensorD **** Dtensors;
         
         //TensorQ's
         TensorQ *** Qtensors;
         
         //Sets everything up for the first solve
         void PreSolve();
         
         //sweepleft
         double sweepleft(const bool change, const int instruction, const bool am_i_master);
         
         //sweepright
         double sweepright(const bool change, const int instruction, const bool am_i_master);
         
         //Load and save functions
         static void MY_HDF5_WRITE(const hid_t file_id, const std::string sPath, Tensor * theTensor);
         static void MY_HDF5_READ( const hid_t file_id, const std::string sPath, Tensor * theTensor);
         void OperatorsOnDisk(const int index, const bool movingRight, const bool store);
         string tempfolder;
         
         void saveMPS(const std::string name, TensorT ** MPSlocation, SyBookkeeper * BKlocation, bool isConverged) const;
         void loadDIM(const std::string name, SyBookkeeper * BKlocation);
         void loadMPS(const std::string name, TensorT ** MPSlocation, bool * isConverged);
         bool makecheckpoints;
         
         //Helper functions for making the boundary operators
         void updateMovingRight(const int index);
         void updateMovingLeft(const int index);
         void deleteTensors(const int index, const bool movingRight);
         void allocateTensors(const int index, const bool movingRight);
         void updateMovingRightSafe(const int cnt);
         void updateMovingRightSafeFirstTime(const int cnt);
         void updateMovingLeftSafe(const int cnt);
         void updateMovingLeftSafe2DM(const int cnt);
         void deleteAllBoundaryOperators();
         static int trianglefunction(const int k, const int glob);
         
         //The storage and functions to handle excited states
         int nStates;
         bool Exc_activated;
         int maxExc;
         double * Exc_Eshifts;
         TensorT *** Exc_MPSs;
         SyBookkeeper ** Exc_BKs;
         TensorO *** Exc_Overlaps;
         double ** prepare_excitations(Sobject * denS);
         void cleanup_excitations(double ** VeffTilde) const;
         void calcVeffTilde(double * result, Sobject * currentS, int state_number);
         void calcOverlapsWithLowerStates();
         
         //For the timings
         double timings[CHEMPS2_TIME_VECLENGTH];
         
   };
}

#endif
