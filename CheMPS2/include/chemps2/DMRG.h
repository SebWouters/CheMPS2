/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2021 Sebastian Wouters

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
#include "TensorOperator.h"
#include "TensorQ.h"
#include "TensorO.h"
#include "Tensor3RDM.h"
#include "TensorGYZ.h"
#include "TensorKM.h"
#include "TwoDM.h"
#include "ThreeDM.h"
#include "Correlations.h"
#include "Heff.h"
#include "Sobject.h"
#include "ConvergenceScheme.h"

#include <hdf5.h>

//For the timings of the different parts of DMRG
#define CHEMPS2_TIME_S_JOIN      0
#define CHEMPS2_TIME_S_SOLVE     1
#define CHEMPS2_TIME_S_SPLIT     2
#define CHEMPS2_TIME_TENS_TOTAL  3
#define CHEMPS2_TIME_TENS_ALLOC  4
#define CHEMPS2_TIME_TENS_FREE   5
#define CHEMPS2_TIME_DISK_WRITE  6
#define CHEMPS2_TIME_DISK_READ   7
#define CHEMPS2_TIME_TENS_CALC   8
#define CHEMPS2_TIME_VECLENGTH   9

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
             \param tmpfolder Temporary folder on a large partition to store the renormalized operators on disk (by default "/tmp")
             \param occupancies ROHF occupancies of a determinant to enlarge in the initial guess in HAM order */
         DMRG(Problem * Probin, ConvergenceScheme * OptSchemeIn, const bool makechkpt=CheMPS2::DMRG_storeMpsOnDisk, const string tmpfolder=CheMPS2::defaultTMPpath, int * occupancies=NULL);
         
         //! Destructor
         virtual ~DMRG();
         
         //! Solver
         /** \return The min. energy encountered so far during the sweeps. */
         double Solve();
         
         //! Reconstruct the renormalized operators when you overwrite the matrix elements with Prob->setMxElement()
         void PreSolve();
         
         //! Calculate the 2-RDM and correlations. Afterwards the MPS is again in LLLLLLLC gauge.
         void calc2DMandCorrelations(){ calc_rdms_and_correlations(false); }
         
         //! Calculate the reduced density matrices and correlations. Afterwards the MPS is again in LLLLLLLC gauge.
         /** \param   do_3rdm Whether or not to calculate the 3-RDM
             \param disk_3rdm Whether or not to use disk in order to avoid storing the full 3-RDM of size L^6 */
         void calc_rdms_and_correlations( const bool do_3rdm, const bool disk_3rdm = false );
         
         //! Get the pointer to the 2-RDM
         /** \return The 2-RDM. Returns a NULL pointer if not yet calculated. */
         TwoDM * get2DM(){ return the2DM; }
         
         //! Get the pointer to the 3-RDM
         /** \return The 3-RDM. Returns a NULL pointer if not yet calculated. */
         ThreeDM * get3DM(){ return the3DM; }

         //! Obtain the symmetrized 4-RDM terms 0.5 * ( Gamma4_ijkl,pqrt + Gamma4_ijkt,pqrl ) with l and t fixed, after the 3-RDM has been calculated.
         /** \param output    Array to store the symmetrized 4-RDM terms in Hamiltonian index notation: output[ i + L * ( j + L * ( k + L * ( p + L * ( q + L * r )))) ] = 0.5 * ( Gamma4_ijkl,pqrt + Gamma4_ijkt,pqrl ).
             \param ham_orb1  The Hamiltonian index of the first  fixed orbital.
             \param ham_orb2  The Hamiltonian index of the second fixed orbital.
             \param last_case If true, everything will be set up to allow to continue sweeping. */
         void Symm4RDM( double * output, const int ham_orb1, const int ham_orb2, const bool last_case );

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
         static void PrintLicense();

         //! Get the number of MPS variables
         int get_num_mps_var() const;
         
      private:
      
         //Setup the DMRG SyBK and MPS (in separate function to allow pushbacks and recreations for excited states)
         void setupBookkeeperAndMPS( int * occupancies=NULL );

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
         
         //The ThreeDM
         ThreeDM * the3DM;
         
         //The Correlations
         Correlations * theCorr;
         
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
         
         //ABCD-tensors
         TensorOperator **** Atensors;
         TensorOperator **** Btensors;
         TensorOperator **** Ctensors;
         TensorOperator **** Dtensors;
         
         //TensorQ's
         TensorQ *** Qtensors;
         
         //Tensors required for the 3-RDM calculation
         Tensor3RDM ***** tensor_3rdm_a_J0_doublet;
         Tensor3RDM ***** tensor_3rdm_a_J1_doublet;
         Tensor3RDM ***** tensor_3rdm_a_J1_quartet;
         Tensor3RDM ***** tensor_3rdm_b_J0_doublet;
         Tensor3RDM ***** tensor_3rdm_b_J1_doublet;
         Tensor3RDM ***** tensor_3rdm_b_J1_quartet;
         Tensor3RDM ***** tensor_3rdm_c_J0_doublet;
         Tensor3RDM ***** tensor_3rdm_c_J1_doublet;
         Tensor3RDM ***** tensor_3rdm_c_J1_quartet;
         Tensor3RDM ***** tensor_3rdm_d_J0_doublet;
         Tensor3RDM ***** tensor_3rdm_d_J1_doublet;
         Tensor3RDM ***** tensor_3rdm_d_J1_quartet;
         
         //Tensors required for the Correlations calculation
         TensorGYZ ** Gtensors;
         TensorGYZ ** Ytensors;
         TensorGYZ ** Ztensors;
         TensorKM  ** Ktensors;
         TensorKM  ** Mtensors;

         // Sweeps
         double sweepleft(  const bool change, const int instruction, const bool am_i_master );
         double sweepright( const bool change, const int instruction, const bool am_i_master );
         double solve_site( const int index, const double dvdson_rtol, const double noise_level, const int virtual_dimension, const bool am_i_master, const bool moving_right, const bool change );

         //Load and save functions
         void MY_HDF5_WRITE_BATCH(const hid_t file_id, const int number, Tensor ** batch, const long long totalsize, const std::string tag);
         void MY_HDF5_READ_BATCH( const hid_t file_id, const int number, Tensor ** batch, const long long totalsize, const std::string tag);
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
         void updateMovingRightSafe2DM(const int cnt);
         void updateMovingLeftSafe(const int cnt);
         void updateMovingLeftSafeFirstTime(const int cnt);
         void updateMovingLeftSafe2DM(const int cnt);
         void deleteAllBoundaryOperators();

         // Helper functions for making the 3-RDM boundary operators
         void update_safe_3rdm_operators( const int boundary );
         void allocate_3rdm_operators( const int boundary );
         void update_3rdm_operators( const int boundary );
         void delete_3rdm_operators( const int boundary );
         
         // Helper functions for making the symmetrized 4-RDM
         void solve_fock( const int dmrg_orb1, const int dmrg_orb2, const double alpha, const double beta );
         static void solve_fock_update_helper( const int index, const int dmrg_orb1, const int dmrg_orb2, const bool moving_right, TensorT ** new_mps, TensorT ** old_mps, SyBookkeeper * new_bk, SyBookkeeper * old_bk, TensorO ** overlaps, TensorL ** regular, TensorL ** trans );
         static void  left_normalize( TensorT * left_mps, TensorT * right_mps );
         static void right_normalize( TensorT * left_mps, TensorT * right_mps );
         void symm_4rdm_helper( double * output, const int ham_orb1, const int ham_orb2, const double alpha, const double beta, const bool add, const double factor );

         //Helper functions for making the Correlations boundary operators
         void update_correlations_tensors(const int siteindex);

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
         void calc_overlaps( const bool moving_right );
         
         // Performance counters
         double timings[ CHEMPS2_TIME_VECLENGTH ];
         long long num_double_write_disk;
         long long num_double_read_disk;
         void print_tensor_update_performance() const;
         
   };
}

#endif
