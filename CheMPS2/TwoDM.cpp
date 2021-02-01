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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <assert.h>

#include <hdf5.h>

#include "TwoDM.h"
#include "Lapack.h"
#include "Options.h"
#include "MPIchemps2.h"
#include "Wigner.h"
#include "Special.h"

using std::max;
using std::cout;
using std::endl;

CheMPS2::TwoDM::TwoDM(const SyBookkeeper * denBKIn, const Problem * ProbIn){

   denBK = denBKIn;
   Prob = ProbIn;
   L = denBK->gL();

   const long long size = ((long long) L ) * ((long long) L ) * ((long long) L ) * ((long long) L );
   assert( INT_MAX >= size );
   two_rdm_A = new double[ size ];
   two_rdm_B = new double[ size ];
   
   //Clear the storage so that an allreduce can be performed in the end
   for (int cnt = 0; cnt < size; cnt++){ two_rdm_A[ cnt ] = 0.0; }
   for (int cnt = 0; cnt < size; cnt++){ two_rdm_B[ cnt ] = 0.0; }

}

CheMPS2::TwoDM::~TwoDM(){

   delete [] two_rdm_A;
   delete [] two_rdm_B;

}

#ifdef CHEMPS2_MPI_COMPILATION
void CheMPS2::TwoDM::mpi_allreduce(){

   const int size = L*L*L*L; // Tested OK in creator TwoDM
   double * temp = new double[ size ];
   MPIchemps2::allreduce_array_double( two_rdm_A, temp, size ); for (int cnt = 0; cnt < size; cnt++){ two_rdm_A[ cnt ] = temp[ cnt ]; }
   MPIchemps2::allreduce_array_double( two_rdm_B, temp, size ); for (int cnt = 0; cnt < size; cnt++){ two_rdm_B[ cnt ] = temp[ cnt ]; }
   delete [] temp;

}
#endif

void CheMPS2::TwoDM::set_2rdm_A_DMRG(const int cnt1, const int cnt2, const int cnt3, const int cnt4, const double value){

   //Prob assumes you use DMRG orbs...
   //Irrep sanity checks are performed in TwoDM::FillSite
   two_rdm_A[ cnt1 + L * ( cnt2 + L * ( cnt3 + L * cnt4 ) ) ] = value;
   two_rdm_A[ cnt2 + L * ( cnt1 + L * ( cnt4 + L * cnt3 ) ) ] = value;
   two_rdm_A[ cnt3 + L * ( cnt4 + L * ( cnt1 + L * cnt2 ) ) ] = value;
   two_rdm_A[ cnt4 + L * ( cnt3 + L * ( cnt2 + L * cnt1 ) ) ] = value;

}

void CheMPS2::TwoDM::set_2rdm_B_DMRG(const int cnt1, const int cnt2, const int cnt3, const int cnt4, const double value){

   //Prob assumes you use DMRG orbs...
   //Irrep sanity checks are performed in TwoDM::FillSite
   two_rdm_B[ cnt1 + L * ( cnt2 + L * ( cnt3 + L * cnt4 ) ) ] = value;
   two_rdm_B[ cnt2 + L * ( cnt1 + L * ( cnt4 + L * cnt3 ) ) ] = value;
   two_rdm_B[ cnt3 + L * ( cnt4 + L * ( cnt1 + L * cnt2 ) ) ] = value;
   two_rdm_B[ cnt4 + L * ( cnt3 + L * ( cnt2 + L * cnt1 ) ) ] = value;

}

double CheMPS2::TwoDM::getTwoDMA_DMRG(const int cnt1, const int cnt2, const int cnt3, const int cnt4) const{

   //Prob assumes you use DMRG orbs...
   const int irrep1 = Prob->gIrrep(cnt1);
   const int irrep2 = Prob->gIrrep(cnt2);
   const int irrep3 = Prob->gIrrep(cnt3);
   const int irrep4 = Prob->gIrrep(cnt4);
   if ( Irreps::directProd(irrep1, irrep2) == Irreps::directProd(irrep3, irrep4) ){
      return two_rdm_A[ cnt1 + L * ( cnt2 + L * ( cnt3 + L * cnt4 ) ) ];
   }
   
   return 0.0;

}

double CheMPS2::TwoDM::getTwoDMB_DMRG(const int cnt1, const int cnt2, const int cnt3, const int cnt4) const{

   //Prob assumes you use DMRG orbs...
   const int irrep1 = Prob->gIrrep(cnt1);
   const int irrep2 = Prob->gIrrep(cnt2);
   const int irrep3 = Prob->gIrrep(cnt3);
   const int irrep4 = Prob->gIrrep(cnt4);
   if ( Irreps::directProd(irrep1, irrep2) == Irreps::directProd(irrep3, irrep4) ){
      return two_rdm_B[ cnt1 + L * ( cnt2 + L * ( cnt3 + L * cnt4 ) ) ];
   }
   
   return 0.0;

}

double CheMPS2::TwoDM::get1RDM_DMRG(const int cnt1, const int cnt2) const{

   //Prob assumes you use DMRG orbs...
   const int irrep1 = Prob->gIrrep(cnt1);
   const int irrep2 = Prob->gIrrep(cnt2);
   if ( irrep1 == irrep2 ){
      double value = 0.0;
      for ( int orbsum = 0; orbsum < L; orbsum++ ){ value += getTwoDMA_DMRG( cnt1, orbsum, cnt2, orbsum ); }
      value = value / ( Prob->gN() - 1.0 );
      return value;
   }
   
   return 0.0;

}

double CheMPS2::TwoDM::spin_density_dmrg( const int cnt1, const int cnt2 ) const{

   //Prob assumes you use DMRG orbs...
   const int irrep1 = Prob->gIrrep(cnt1);
   const int irrep2 = Prob->gIrrep(cnt2);
   if ( irrep1 == irrep2 ){
      const int two_s = Prob->gTwoS();
      if ( two_s > 0 ){
         double value = ( 2 - Prob->gN() ) * get1RDM_DMRG( cnt1, cnt2 );
         for ( int orb = 0; orb < Prob->gL(); orb++ ){
            value -= ( getTwoDMA_DMRG( cnt1, orb, orb, cnt2 ) + getTwoDMB_DMRG( cnt1, orb, orb, cnt2 ) );
         }
         value = 1.5 * value / ( 0.5 * two_s + 1 );
         return value;
      }
   }

   return 0.0;

}

double CheMPS2::TwoDM::getTwoDMA_HAM(const int cnt1, const int cnt2, const int cnt3, const int cnt4) const{

   //Prob assumes you use DMRG orbs... f1 converts HAM orbs to DMRG orbs
   if ( Prob->gReorder() ){
      return getTwoDMA_DMRG( Prob->gf1(cnt1), Prob->gf1(cnt2), Prob->gf1(cnt3), Prob->gf1(cnt4) );
   }
   return getTwoDMA_DMRG( cnt1, cnt2, cnt3, cnt4 );

}

double CheMPS2::TwoDM::getTwoDMB_HAM(const int cnt1, const int cnt2, const int cnt3, const int cnt4) const{

   //Prob assumes you use DMRG orbs... f1 converts HAM orbs to DMRG orbs
   if ( Prob->gReorder() ){
      return getTwoDMB_DMRG( Prob->gf1(cnt1), Prob->gf1(cnt2), Prob->gf1(cnt3), Prob->gf1(cnt4) );
   }
   return getTwoDMB_DMRG( cnt1, cnt2, cnt3, cnt4 );

}

double CheMPS2::TwoDM::get1RDM_HAM(const int cnt1, const int cnt2) const{

   //Prob assumes you use DMRG orbs... f1 converts HAM orbs to DMRG orbs
   if ( Prob->gReorder() ){
      return get1RDM_DMRG( Prob->gf1(cnt1), Prob->gf1(cnt2) );
   }
   return get1RDM_DMRG( cnt1, cnt2 );

}

double CheMPS2::TwoDM::spin_density_ham( const int cnt1, const int cnt2 ) const{

   //Prob assumes you use DMRG orbs... f1 converts HAM orbs to DMRG orbs
   if ( Prob->gReorder() ){
      return spin_density_dmrg( Prob->gf1(cnt1), Prob->gf1(cnt2) );
   }
   return spin_density_dmrg( cnt1, cnt2 );

}

double CheMPS2::TwoDM::trace() const{

   double val = 0.0;
   for (int cnt1=0; cnt1<L; cnt1++){
      for (int cnt2=0; cnt2<L; cnt2++){
         val += getTwoDMA_DMRG(cnt1,cnt2,cnt1,cnt2);
      }
   }
   return val;

}

double CheMPS2::TwoDM::energy() const{

   double val = 0.0;
   for (int cnt1=0; cnt1<L; cnt1++){
      for (int cnt2=0; cnt2<L; cnt2++){
         for (int cnt3=0; cnt3<L; cnt3++){
            for (int cnt4=0; cnt4<L; cnt4++){
               val += getTwoDMA_DMRG(cnt1,cnt2,cnt3,cnt4) * Prob->gMxElement(cnt1,cnt2,cnt3,cnt4);
            }
         }
      }
   }
   val *= 0.5;
   return val + Prob->gEconst();

}

void CheMPS2::TwoDM::print_noon() const{

   int lwork = 3 * L;
   double * OneRDM = new double[ L * L ];
   double * work   = new double[ lwork ];
   double * eigs   = new double[ L ];

   Irreps my_irreps( Prob->gSy() );

   for ( int irrep = 0; irrep < denBK->getNumberOfIrreps(); irrep++ ){
   
      int jump1 = 0;
      for ( int orb1 = 0; orb1 < L; orb1++ ){
         if ( Prob->gIrrep( orb1 ) == irrep ){
            int jump2 = jump1;
            for ( int orb2 = orb1; orb2 < L; orb2++ ){
               if ( Prob->gIrrep( orb2 ) == irrep ){
                  const double value = get1RDM_DMRG( orb1, orb2 );
                  OneRDM[ jump1 + L * jump2 ] = value;
                  OneRDM[ jump2 + L * jump1 ] = value;
                  jump2 += 1;
               }
            }
            jump1 += 1;
         }
      }
      
      if ( jump1 > 0 ){
         char jobz = 'N'; // Eigenvalues only
         char uplo = 'U';
         int lda = L;
         int info;
         dsyev_(&jobz, &uplo, &jump1, OneRDM, &lda, eigs, work, &lwork, &info);
         cout << "   NOON of irrep " << my_irreps.getIrrepName( irrep ) << " = [ ";
         for ( int cnt = 0; cnt < jump1 - 1; cnt++ ){ cout << eigs[ jump1 - 1 - cnt ] << " , "; } // Print from large to small
         cout << eigs[ 0 ] << " ]." << endl;
      }
   
   }
   delete [] OneRDM;
   delete [] work;
   delete [] eigs;

}

void CheMPS2::TwoDM::save_HAM( const string filename ) const{

   // Create an array with the 2-RDM in the ORIGINAL HAM indices
   const int total_size = L * L * L * L;
   double * local_array = new double[ total_size ];
   for ( int ham4 = 0; ham4 < L; ham4++ ){
      for ( int ham3 = 0; ham3 < L; ham3++ ){
         for ( int ham2 = 0; ham2 < L; ham2++ ){
            for ( int ham1 = 0; ham1 < L; ham1++ ){
                local_array[ ham1 + L * ( ham2 + L * ( ham3 + L * ham4 ) ) ] = getTwoDMA_HAM( ham1, ham2, ham3, ham4 );
            }
         }
      }
   }

   // (Re)create the HDF5 file with the 2-RDM
   hid_t   file_id      = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
   hsize_t dimarray     = total_size;
   hid_t   group_id     = H5Gcreate( file_id, "2-RDM", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
   hid_t   dataspace_id = H5Screate_simple( 1, &dimarray, NULL );
   hid_t   dataset_id   = H5Dcreate( group_id, "elements", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

   H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, local_array );

   H5Dclose( dataset_id );
   H5Sclose( dataspace_id );
   H5Gclose( group_id );
   H5Fclose( file_id );

   // Deallocate the array
   delete [] local_array;

   cout << "Saved the 2-RDM to the file " << filename << endl;

}

void CheMPS2::TwoDM::save() const{

   hid_t   file_id  = H5Fcreate(CheMPS2::TWO_RDM_storagename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   hsize_t dimarray = L*L*L*L;
   {
      hid_t group_id = H5Gcreate(file_id, "two_rdm_A", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

         hid_t dataspace_id     = H5Screate_simple(1, &dimarray, NULL);
         hid_t dataset_id       = H5Dcreate(group_id, "elements", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, two_rdm_A);

         H5Dclose(dataset_id);
         H5Sclose(dataspace_id);

      H5Gclose(group_id);
   }
   {
      hid_t group_id = H5Gcreate(file_id, "two_rdm_B", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

         hid_t dataspace_id     = H5Screate_simple(1, &dimarray, NULL);
         hid_t dataset_id       = H5Dcreate(group_id, "elements", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, two_rdm_B);
         H5Dclose(dataset_id);
         H5Sclose(dataspace_id);

      H5Gclose(group_id);
   }
   H5Fclose(file_id);

}

void CheMPS2::TwoDM::read(){

   hid_t file_id = H5Fopen(CheMPS2::TWO_RDM_storagename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   {
      hid_t group_id = H5Gopen(file_id, "two_rdm_A", H5P_DEFAULT);

         hid_t dataset_id = H5Dopen(group_id, "elements", H5P_DEFAULT);
         H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, two_rdm_A);
         H5Dclose(dataset_id);

      H5Gclose(group_id);
   }
   {
      hid_t group_id = H5Gopen(file_id, "two_rdm_B", H5P_DEFAULT);

         hid_t dataset_id = H5Dopen(group_id, "elements", H5P_DEFAULT);
         H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, two_rdm_B);
         H5Dclose(dataset_id);

      H5Gclose(group_id);
   }
   H5Fclose(file_id);
   
   std::cout << "TwoDM::read : Everything loaded!" << std::endl;

}

void CheMPS2::TwoDM::write2DMAfile(const string filename) const{

   int * psi2molpro = new int[ denBK->getNumberOfIrreps() ];
   Irreps my_irreps( Prob->gSy() );
   my_irreps.symm_psi2molpro( psi2molpro );
   
   FILE * capturing;
   capturing = fopen( filename.c_str(), "w" ); // "w" with fopen means truncate file
   fprintf( capturing, " &2-RDM NORB= %d,NELEC= %d,MS2= %d,\n", L, Prob->gN(), Prob->gTwoS() );
   fprintf( capturing, "  ORBSYM=" );
   for (int HamOrb=0; HamOrb<L; HamOrb++){
      const int DMRGOrb = ( Prob->gReorder() ) ? Prob->gf1( HamOrb ) : HamOrb;
      fprintf( capturing, "%d,", psi2molpro[ Prob->gIrrep( DMRGOrb ) ] );
   }
   fprintf( capturing, "\n  ISYM=%d,\n /\n", psi2molpro[ Prob->gIrrep() ] );
   delete [] psi2molpro;
   
   for (int ham_p=0; ham_p<L; ham_p++){
      const int dmrg_p = (Prob->gReorder())?Prob->gf1(ham_p):ham_p;
      for (int ham_q=0; ham_q<=ham_p; ham_q++){ // p >= q
         const int dmrg_q = (Prob->gReorder())?Prob->gf1(ham_q):ham_q;
         const int irrep_pq = Irreps::directProd( Prob->gIrrep(dmrg_p), Prob->gIrrep(dmrg_q) );
         for (int ham_r=0; ham_r<=ham_p; ham_r++){ // p >= r
            const int dmrg_r = (Prob->gReorder())?Prob->gf1(ham_r):ham_r;
            for (int ham_s=0; ham_s<=ham_p; ham_s++){ // p >= s
               const int dmrg_s = (Prob->gReorder())?Prob->gf1(ham_s):ham_s;
               const int irrep_rs = Irreps::directProd( Prob->gIrrep(dmrg_r), Prob->gIrrep(dmrg_s) );
               if ( irrep_pq == irrep_rs ){
                  const int num_equal = (( ham_q == ham_p ) ? 1 : 0 )
                                      + (( ham_r == ham_p ) ? 1 : 0 )
                                      + (( ham_s == ham_p ) ? 1 : 0 );
                  /*   1. p > q,r,s
                       2. p==q > r,s
                       3. p==r > q,s
                       4. p==s > q,r
                       5. p==q==r > s
                       6. p==q==s > r
                       7. p==r==s > q
                       8. p==r==s==q                  
                  While 2-4 are inequivalent ( num_equal == 1 ), 5-7 are equivalent ( num_equal == 2 ). Hence:  */
                  if ( ( num_equal != 2 ) || ( ham_p > ham_s ) ){
                     const double value = getTwoDMA_DMRG(dmrg_p, dmrg_r, dmrg_q, dmrg_s);
                     fprintf( capturing, " % 23.16E %3d %3d %3d %3d\n", value, ham_p+1, ham_q+1, ham_r+1, ham_s+1 );
                  }
               }
            }
         }
      }
   }
   
   // 1-RDM in Hamiltonian indices with p >= q
   const double prefactor = 1.0 / ( Prob->gN() - 1.0 );
   for (int ham_p=0; ham_p<L; ham_p++){
      const int dmrg_p = (Prob->gReorder())?Prob->gf1(ham_p):ham_p;
      for (int ham_q=0; ham_q<=ham_p; ham_q++){
         const int dmrg_q = (Prob->gReorder())?Prob->gf1(ham_q):ham_q;
         if ( Prob->gIrrep(dmrg_p) == Prob->gIrrep(dmrg_q) ){
            double value = 0.0;
            for ( int orbsum = 0; orbsum < L; orbsum++ ){ value += getTwoDMA_DMRG( dmrg_p, orbsum, dmrg_q, orbsum ); }
            value *= prefactor;
            fprintf( capturing, " % 23.16E %3d %3d %3d %3d\n", value, ham_p+1, ham_q+1, 0, 0 );
         }
      }
   }
   
   // 0-RDM == Norm of the wavefunction?
   fprintf( capturing, " % 23.16E %3d %3d %3d %3d", 1.0, 0, 0, 0, 0 );
   fclose( capturing );
   cout << "Created the file " << filename << "." << endl;

}

void CheMPS2::TwoDM::FillSite(TensorT * denT, TensorL *** Ltens, TensorF0 **** F0tens, TensorF1 **** F1tens, TensorS0 **** S0tens, TensorS1 **** S1tens){

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   const int theindex = denT->gIndex();
   const int DIM = max(denBK->gMaxDimAtBound(theindex), denBK->gMaxDimAtBound(theindex+1));
   
   #ifdef CHEMPS2_MPI_COMPILATION
   if ( MPIRANK == MPI_CHEMPS2_MASTER )
   #endif
   {
      //Diagram 1
      const double d1 = doD1(denT);
      set_2rdm_A_DMRG(theindex,theindex,theindex,theindex, 2*d1);
      set_2rdm_B_DMRG(theindex,theindex,theindex,theindex,-2*d1);
   }
   
   #pragma omp parallel
   {
   
      double * workmem  = new double[DIM*DIM];
      double * workmem2 = new double[DIM*DIM];
      
      #pragma omp for schedule(static) nowait
      for (int j_index=theindex+1; j_index<L; j_index++){
         if (denBK->gIrrep(j_index) == denBK->gIrrep(theindex)){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_q( L, j_index ) ) //Everyone owns the L-tensors --> task division based on Q-tensor ownership
            #endif
            {
               //Diagram 2
               const double d2 = doD2(denT, Ltens[theindex][j_index-theindex-1], workmem);
               set_2rdm_A_DMRG(theindex,j_index,theindex,theindex, 2*d2);
               set_2rdm_B_DMRG(theindex,j_index,theindex,theindex,-2*d2);
            }
         }
      }
      
      const int dimTriangle = L - theindex - 1;
      const int upperboundTriangle = ( dimTriangle * ( dimTriangle + 1 ) ) / 2;
      int result[ 2 ];
      #pragma omp for schedule(static) nowait
      for ( int global = 0; global < upperboundTriangle; global++ ){
         Special::invert_triangle_two( global, result );
         const int j_index = L - 1 - result[ 1 ];
         const int k_index = j_index + result[ 0 ];
         if ( denBK->gIrrep( j_index ) == denBK->gIrrep( k_index )){

            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_absigma( j_index, k_index ) )
            #endif
            {
               //Diagram 3
               const double d3 = doD3(denT, S0tens[theindex][k_index-j_index][j_index-theindex-1], workmem);
               set_2rdm_A_DMRG(theindex,theindex,j_index,k_index, 2*d3);
               set_2rdm_B_DMRG(theindex,theindex,j_index,k_index,-2*d3);
            }

            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_cdf( L, j_index, k_index ) )
            #endif
            {
               //Diagrams 4,5 & 6
               const double d4 = doD4(denT, F0tens[theindex][k_index-j_index][j_index-theindex-1], workmem);
               const double d5 = doD5(denT, F0tens[theindex][k_index-j_index][j_index-theindex-1], workmem);
               const double d6 = doD6(denT, F1tens[theindex][k_index-j_index][j_index-theindex-1], workmem);
               set_2rdm_A_DMRG(theindex,j_index,k_index,theindex, -2*d4 - 2*d5 - 3*d6);
               set_2rdm_B_DMRG(theindex,j_index,k_index,theindex, -2*d4 - 2*d5 +   d6);
               set_2rdm_A_DMRG(theindex,j_index,theindex,k_index,  4*d4 + 4*d5);
               set_2rdm_B_DMRG(theindex,j_index,theindex,k_index,  2*d6);
            }
         }
      }

      #pragma omp for schedule(static) nowait
      for (int g_index=0; g_index<theindex; g_index++){
         if (denBK->gIrrep(g_index) == denBK->gIrrep(theindex)){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_q( L, g_index ) ) //Everyone owns the L-tensors --> task division based on Q-tensor ownership
            #endif
            {
               //Diagram 7
               const double d7 = doD7(denT, Ltens[theindex-1][theindex-g_index-1], workmem);
               set_2rdm_A_DMRG(g_index,theindex,theindex,theindex, 2*d7);
               set_2rdm_B_DMRG(g_index,theindex,theindex,theindex,-2*d7);
            }
         }
      }

      const int globalsize8to12 = theindex * ( L - 1 - theindex );
      #pragma omp for schedule(static) nowait
      for (int gj_index=0; gj_index<globalsize8to12; gj_index++){
         const int g_index = gj_index % theindex;
         const int j_index = ( gj_index / theindex ) + theindex + 1;
         const int I_g = denBK->gIrrep(g_index);
         if (denBK->gIrrep(g_index) == denBK->gIrrep(j_index)){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_absigma( g_index, j_index ) ) //Everyone owns the L-tensors --> task division based on ABSigma-tensor ownership
            #endif
            {
               //Diagrams 8,9,10 & 11
               const double d8 = doD8(denT, Ltens[theindex-1][theindex-g_index-1], Ltens[theindex][j_index-theindex-1], workmem, workmem2, I_g);
               double d9, d10, d11;
               doD9andD10andD11(denT, Ltens[theindex-1][theindex-g_index-1], Ltens[theindex][j_index-theindex-1], workmem, workmem2, &d9, &d10, &d11, I_g);
               set_2rdm_A_DMRG(g_index,theindex,j_index,theindex, -4*d8-d9);
               set_2rdm_A_DMRG(g_index,theindex,theindex,j_index, 2*d8 + d11);
               set_2rdm_B_DMRG(g_index,theindex,j_index,theindex, d9 - 2*d10);
               set_2rdm_B_DMRG(g_index,theindex,theindex,j_index, 2*d8 + 2*d10 - d11);
               
               //Diagram 12
               const double d12 = doD12(denT, Ltens[theindex-1][theindex-g_index-1], Ltens[theindex][j_index-theindex-1], workmem, workmem2, I_g);
               set_2rdm_A_DMRG(g_index,j_index,theindex,theindex, 2*d12);
               set_2rdm_B_DMRG(g_index,j_index,theindex,theindex,-2*d12);
            }
         }
      }

      const int globalsize = theindex * upperboundTriangle;
      #pragma omp for schedule(static) nowait
      for ( int gjk_index = 0; gjk_index < globalsize; gjk_index++ ){
         const int g_index = gjk_index % theindex;
         const int global  = gjk_index / theindex;
         Special::invert_triangle_two( global, result );
         const int j_index = L - 1 - result[ 1 ];
         const int k_index = j_index + result[ 0 ];
         const int I_g = denBK->gIrrep( g_index );
         const int cnt1 = k_index - j_index;
         const int cnt2 = j_index - theindex - 1;

         if (Irreps::directProd(I_g, denBK->gIrrep(theindex)) == Irreps::directProd(denBK->gIrrep(j_index), denBK->gIrrep(k_index))){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_absigma( j_index, k_index ) )
            #endif
            {
               //Diagrams 13,14,15 & 16
               const double d13 = doD13(denT, Ltens[theindex-1][theindex-g_index-1], S0tens[theindex][cnt1][cnt2], workmem, workmem2, I_g);
               const double d14 = doD14(denT, Ltens[theindex-1][theindex-g_index-1], S0tens[theindex][cnt1][cnt2], workmem, workmem2, I_g);
               double d15 = 0.0;
               double d16 = 0.0;
               if (k_index>j_index){
                  d15 = doD15(denT, Ltens[theindex-1][theindex-g_index-1], S1tens[theindex][cnt1][cnt2], workmem, workmem2, I_g);
                  d16 = doD16(denT, Ltens[theindex-1][theindex-g_index-1], S1tens[theindex][cnt1][cnt2], workmem, workmem2, I_g);
               }
               set_2rdm_A_DMRG(g_index,theindex,j_index,k_index, 2*d13 + 2*d14 + 3*d15 + 3*d16);
               set_2rdm_A_DMRG(g_index,theindex,k_index,j_index, 2*d13 + 2*d14 - 3*d15 - 3*d16);
               set_2rdm_B_DMRG(g_index,theindex,j_index,k_index,-2*d13 - 2*d14 +   d15 +   d16);
               set_2rdm_B_DMRG(g_index,theindex,k_index,j_index,-2*d13 - 2*d14 -   d15 -   d16);
            }
            
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_cdf( L, j_index, k_index ) )
            #endif
            {
               //Diagrams 17,18,19 & 20
               const double d17 = doD17orD21(denT, Ltens[theindex-1][theindex-g_index-1], F0tens[theindex][cnt1][cnt2], workmem, workmem2, I_g, true);
               const double d18 = doD18orD22(denT, Ltens[theindex-1][theindex-g_index-1], F0tens[theindex][cnt1][cnt2], workmem, workmem2, I_g, true);
               const double d19 = doD19orD23(denT, Ltens[theindex-1][theindex-g_index-1], F1tens[theindex][cnt1][cnt2], workmem, workmem2, I_g, true);
               const double d20 = doD20orD24(denT, Ltens[theindex-1][theindex-g_index-1], F1tens[theindex][cnt1][cnt2], workmem, workmem2, I_g, true);
               set_2rdm_A_DMRG(g_index,j_index,k_index,theindex, -2*d17 - 2*d18 - 3*d19 - 3*d20);
               set_2rdm_A_DMRG(g_index,j_index,theindex,k_index,  4*d17 + 4*d18                );
               set_2rdm_B_DMRG(g_index,j_index,k_index,theindex, -2*d17 - 2*d18 +   d19 +   d20);
               set_2rdm_B_DMRG(g_index,j_index,theindex,k_index,                  2*d19 + 2*d20);
               
               //Diagrams 21,22,23 & 24
               const double d21 = doD17orD21(denT, Ltens[theindex-1][theindex-g_index-1], F0tens[theindex][cnt1][cnt2], workmem, workmem2, I_g, false);
               const double d22 = doD18orD22(denT, Ltens[theindex-1][theindex-g_index-1], F0tens[theindex][cnt1][cnt2], workmem, workmem2, I_g, false);
               const double d23 = doD19orD23(denT, Ltens[theindex-1][theindex-g_index-1], F1tens[theindex][cnt1][cnt2], workmem, workmem2, I_g, false);
               const double d24 = doD20orD24(denT, Ltens[theindex-1][theindex-g_index-1], F1tens[theindex][cnt1][cnt2], workmem, workmem2, I_g, false);
               set_2rdm_A_DMRG(g_index,k_index,j_index,theindex, -2*d21 - 2*d22 - 3*d23 - 3*d24);
               set_2rdm_A_DMRG(g_index,k_index,theindex,j_index,  4*d21 + 4*d22                );
               set_2rdm_B_DMRG(g_index,k_index,j_index,theindex, -2*d21 - 2*d22 +   d23 +   d24);
               set_2rdm_B_DMRG(g_index,k_index,theindex,j_index,                  2*d23 + 2*d24);
            }
         }
      }

      delete [] workmem;
      delete [] workmem2;
   
   }

}

void CheMPS2::TwoDM::correct_higher_multiplicities(){

   if ( Prob->gTwoS() != 0 ){
      double alpha = 1.0 / ( Prob->gTwoS() + 1.0 );
      int length   = L*L*L*L;
      int inc      = 1;
      dscal_( &length, &alpha, two_rdm_A, &inc );
      dscal_( &length, &alpha, two_rdm_B, &inc );
   }

}

double CheMPS2::TwoDM::doD1(TensorT * denT){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSL = denBK->gTwoSmin(theindex,NL); TwoSL<= denBK->gTwoSmax(theindex,NL); TwoSL+=2){
         for (int IL = 0; IL<denBK->getNumberOfIrreps(); IL++){
            int dimL = denBK->gCurrentDim(theindex,  NL,  TwoSL,IL);
            int dimR = denBK->gCurrentDim(theindex+1,NL+2,TwoSL,IL);
            if ((dimL>0) && (dimR>0)){
            
               double * Tblock = denT->gStorage(NL,TwoSL,IL,NL+2,TwoSL,IL);
               
               int length = dimL*dimR;
               int inc = 1;
               total += (TwoSL+1) * ddot_(&length, Tblock, &inc, Tblock, &inc);
               
            }
         }
      }
   }
   
   return total;

}

double CheMPS2::TwoDM::doD2(TensorT * denT, TensorL * Lright, double * workmem){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSL = denBK->gTwoSmin(theindex,NL); TwoSL<= denBK->gTwoSmax(theindex,NL); TwoSL+=2){
         for (int IL = 0; IL<denBK->getNumberOfIrreps(); IL++){
            for (int TwoSR = TwoSL-1; TwoSR<=TwoSL+1; TwoSR+=2){
            
               int IRup = Irreps::directProd(IL, denBK->gIrrep(theindex)); 
               
               int dimL     = denBK->gCurrentDim(theindex,  NL,  TwoSL,IL);
               int dimRdown = denBK->gCurrentDim(theindex+1,NL+2,TwoSL,IL);
               int dimRup   = denBK->gCurrentDim(theindex+1,NL+1,TwoSR,IRup);
            
               if ((dimL>0) && (dimRup>0) && (dimRdown>0)){
               
                  double * Tdown = denT->gStorage(NL,TwoSL,IL,NL+2,TwoSL,IL  );
                  double * Tup   = denT->gStorage(NL,TwoSL,IL,NL+1,TwoSR,IRup);
                  double * Lblock = Lright->gStorage(NL+1,TwoSR,IRup,NL+2,TwoSL,IL);
               
                  char trans = 'T';
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0; //set
                  dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,Tdown,&dimL,Lblock,&dimRup,&beta,workmem,&dimL);

                  const double factor = Special::phase( TwoSL + 1 - TwoSR ) * 0.5 * sqrt((TwoSL+1)*(TwoSR+1.0));

                  int length = dimL * dimRup;
                  int inc = 1;
                  total += factor * ddot_(&length, workmem, &inc, Tup, &inc);

               }
            }
         }
      }
   }
   
   return total;

}

double CheMPS2::TwoDM::doD3(TensorT * denT, TensorS0 * S0right, double * workmem){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSL = denBK->gTwoSmin(theindex,NL); TwoSL<= denBK->gTwoSmax(theindex,NL); TwoSL+=2){
         for (int IL = 0; IL<denBK->getNumberOfIrreps(); IL++){
               
            int dimL     = denBK->gCurrentDim(theindex,  NL,  TwoSL,IL);
            int dimRdown = denBK->gCurrentDim(theindex+1,NL,  TwoSL,IL);
            int dimRup   = denBK->gCurrentDim(theindex+1,NL+2,TwoSL,IL);
            
            if ((dimL>0) && (dimRup>0) && (dimRdown>0)){
               
               double * Tdown = denT->gStorage(NL,TwoSL,IL,NL,  TwoSL,IL);
               double * Tup   = denT->gStorage(NL,TwoSL,IL,NL+2,TwoSL,IL);
               double * S0block = S0right->gStorage(NL, TwoSL, IL, NL+2, TwoSL, IL);
               
               char notrans = 'N';
               double alpha = 1.0;
               double beta = 0.0; //set
               dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,Tdown,&dimL,S0block,&dimRdown,&beta,workmem,&dimL);
               
               double factor = sqrt(0.5) * (TwoSL+1);
               
               int length = dimL * dimRup;
               int inc = 1;
               total += factor * ddot_(&length, workmem, &inc, Tup, &inc);

            }
         }
      }
   }

   return total;

}

double CheMPS2::TwoDM::doD4(TensorT * denT, TensorF0 * F0right, double * workmem){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSL = denBK->gTwoSmin(theindex,NL); TwoSL<= denBK->gTwoSmax(theindex,NL); TwoSL+=2){
         for (int IL = 0; IL<denBK->getNumberOfIrreps(); IL++){
               
            int dimL     = denBK->gCurrentDim(theindex,  NL,  TwoSL,IL);
            int dimR     = denBK->gCurrentDim(theindex+1,NL+2,TwoSL,IL);
            
            if ((dimL>0) && (dimR>0)){
               
               double * Tblock  =    denT->gStorage(NL,   TwoSL, IL, NL+2, TwoSL, IL);
               double * F0block = F0right->gStorage(NL+2, TwoSL, IL, NL+2, TwoSL, IL);
               
               char notrans = 'N';
               double alpha = 1.0;
               double beta = 0.0; //set
               dgemm_(&notrans,&notrans,&dimL,&dimR,&dimR,&alpha,Tblock,&dimL,F0block,&dimR,&beta,workmem,&dimL);
               
               double factor = sqrt(0.5) * (TwoSL+1);
               
               int length = dimL * dimR;
               int inc = 1;
               total += factor * ddot_(&length, workmem, &inc, Tblock, &inc);

            }
         }
      }
   }

   return total;

}

double CheMPS2::TwoDM::doD5(TensorT * denT, TensorF0 * F0right, double * workmem){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSL = denBK->gTwoSmin(theindex,NL); TwoSL<= denBK->gTwoSmax(theindex,NL); TwoSL+=2){
         for (int IL = 0; IL<denBK->getNumberOfIrreps(); IL++){
            for (int TwoSR=TwoSL-1; TwoSR<=TwoSL+1; TwoSR+=2){
               
               int IR = Irreps::directProd(IL,denBK->gIrrep(theindex));
               int dimL     = denBK->gCurrentDim(theindex,  NL,  TwoSL,IL);
               int dimR     = denBK->gCurrentDim(theindex+1,NL+1,TwoSR,IR);
            
               if ((dimL>0) && (dimR>0)){
               
                  double * Tblock  =    denT->gStorage(NL,   TwoSL, IL, NL+1, TwoSR, IR);
                  double * F0block = F0right->gStorage(NL+1, TwoSR, IR, NL+1, TwoSR, IR);
               
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0; //set
                  dgemm_(&notrans,&notrans,&dimL,&dimR,&dimR,&alpha,Tblock,&dimL,F0block,&dimR,&beta,workmem,&dimL);
               
                  double factor = 0.5 * sqrt(0.5) * (TwoSR+1);
                  
                  int length = dimL * dimR;
                  int inc = 1;
                  total += factor * ddot_(&length, workmem, &inc, Tblock, &inc);
                  
               }
            }
         }
      }
   }

   return total;

}

double CheMPS2::TwoDM::doD6(TensorT * denT, TensorF1 * F1right, double * workmem){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSL = denBK->gTwoSmin(theindex,NL); TwoSL<= denBK->gTwoSmax(theindex,NL); TwoSL+=2){
         for (int IL = 0; IL<denBK->getNumberOfIrreps(); IL++){
            for (int TwoSRup=TwoSL-1; TwoSRup<=TwoSL+1; TwoSRup+=2){
               for (int TwoSRdown=TwoSL-1; TwoSRdown<=TwoSL+1; TwoSRdown+=2){
               
                  int IR = Irreps::directProd(IL,denBK->gIrrep(theindex));
                  int dimL     = denBK->gCurrentDim(theindex,  NL,  TwoSL,    IL);
                  int dimRup   = denBK->gCurrentDim(theindex+1,NL+1,TwoSRup,  IR);
                  int dimRdown = denBK->gCurrentDim(theindex+1,NL+1,TwoSRdown,IR);
            
                  if ((dimL>0) && (dimRup>0) && (dimRdown>0)){
               
                     double * Tup   =      denT->gStorage(NL,   TwoSL,     IL, NL+1, TwoSRup,   IR);
                     double * Tdown =      denT->gStorage(NL,   TwoSL,     IL, NL+1, TwoSRdown, IR);
                     double * F1block = F1right->gStorage(NL+1, TwoSRdown, IR, NL+1, TwoSRup,   IR);
               
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta = 0.0; //set
                     dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,Tdown,&dimL,F1block,&dimRdown,&beta,workmem,&dimL);

                     const double factor = sqrt((TwoSRup+1)/3.0) * ( TwoSRdown + 1 )
                                         * Special::phase( TwoSL + TwoSRdown - 1 )
                                         * Wigner::wigner6j( 1, 1, 2, TwoSRup, TwoSRdown, TwoSL );

                     int length = dimL * dimRup;
                     int inc = 1;
                     total += factor * ddot_(&length, workmem, &inc, Tup, &inc);
                     
                  }
               }
            }
         }
      }
   }

   return total;

}

double CheMPS2::TwoDM::doD7(TensorT * denT, TensorL * Lleft, double * workmem){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NLup = denBK->gNmin(theindex); NLup<=denBK->gNmax(theindex); NLup++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NLup); TwoSLup<= denBK->gTwoSmax(theindex,NLup); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex, NLup, TwoSLup, ILup);
         
            for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
               
               int IR = Irreps::directProd(ILup,denBK->gIrrep(theindex));
               int dimLdown = denBK->gCurrentDim(theindex,   NLup-1, TwoSLdown, IR);
               int dimR     = denBK->gCurrentDim(theindex+1, NLup+1, TwoSLdown, IR);
            
               if ((dimLup>0) && (dimLdown>0) && (dimR>0)){
               
                  double * Tup   =   denT->gStorage(NLup,   TwoSLup,   ILup, NLup+1, TwoSLdown, IR);
                  double * Tdown =   denT->gStorage(NLup-1, TwoSLdown, IR,   NLup+1, TwoSLdown, IR);
                  double * Lblock = Lleft->gStorage(NLup-1, TwoSLdown, IR,   NLup,   TwoSLup,   ILup);
               
                  char trans = 'T';
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0; //set
                  dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&alpha,Lblock,&dimLdown,Tdown,&dimLdown,&beta,workmem,&dimLup);

                  const double factor = 0.5 * sqrt((TwoSLdown+1)*(TwoSLup+1.0)) * Special::phase( TwoSLup - TwoSLdown + 3 );

                  int length = dimLup * dimR;
                  int inc = 1;
                  total += factor * ddot_(&length, workmem, &inc, Tup, &inc);

               }
            }
         }
      }
   }

   return total;

}

double CheMPS2::TwoDM::doD8(TensorT * denT, TensorL * Lleft, TensorL * Lright, double * workmem, double * workmem2, int Irrep_g){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            int dimRup = denBK->gCurrentDim(theindex+1, NL+2, TwoSLup, ILup);
            
            if ((dimLup>0) && (dimRup>0)){
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  
                  int Idown = Irreps::directProd(ILup,Irrep_g);
                  
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, Idown);
                  int dimRdown = denBK->gCurrentDim(theindex+1, NL+1, TwoSLdown, Idown);
               
                  if ((dimLdown>0) && (dimRdown>0)){
                  
                     double * Tup       =   denT->gStorage(NL,   TwoSLup,   ILup,  NL+2, TwoSLup,   ILup);
                     double * Tdown     =   denT->gStorage(NL-1, TwoSLdown, Idown, NL+1, TwoSLdown, Idown);
                     double * LleftBlk  =  Lleft->gStorage(NL-1, TwoSLdown, Idown, NL,   TwoSLup, ILup);
                     double * LrightBlk = Lright->gStorage(NL+1, TwoSLdown, Idown, NL+2, TwoSLup, ILup);
                  
                     char trans = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta = 0.0; //set
                     dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,LleftBlk,&dimLdown,Tdown,&dimLdown,&beta,workmem,&dimLup);
                  
                     dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,LrightBlk,&dimRdown,&beta,workmem2,&dimLup);
                  
                     double factor = -0.5 * (TwoSLup+1);
                     
                     int length = dimLup * dimRup;
                     int inc = 1;
                     total += factor * ddot_(&length, workmem2, &inc, Tup, &inc);
                  
                  }
               }
            }
         }
      }
   }

   return total;

}

void CheMPS2::TwoDM::doD9andD10andD11(TensorT * denT, TensorL * Lleft, TensorL * Lright, double * workmem, double * workmem2, double * d9, double * d10, double * d11, int Irrep_g){

   d9[0]  = 0.0;
   d10[0] = 0.0;
   d11[0] = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex, NL, TwoSLup, ILup);
            if (dimLup>0){
            
               int IRup   = Irreps::directProd(ILup,   denBK->gIrrep(theindex));
               int ILdown = Irreps::directProd(ILup,   Irrep_g);
               int IRdown = Irreps::directProd(ILdown, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  for (int TwoSRup=TwoSLup-1; TwoSRup<=TwoSLup+1; TwoSRup+=2){
                     for (int TwoSRdown=TwoSRup-1; TwoSRdown<=TwoSRup+1; TwoSRdown+=2){
                        if ((TwoSLdown>=0) && (TwoSRup>=0) && (TwoSRdown>=0) && (abs(TwoSLdown - TwoSRdown)<=1)){
                        
                           int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                           int dimRup   = denBK->gCurrentDim(theindex+1, NL+1, TwoSRup,   IRup);
                           int dimRdown = denBK->gCurrentDim(theindex+1, NL,   TwoSRdown, IRdown);
               
                           if ((dimLdown>0) && (dimRup>0) && (dimRdown>0)){
                           
                              double * T_up      =   denT->gStorage(NL,   TwoSLup,   ILup,   NL+1, TwoSRup,   IRup);
                              double * T_down    =   denT->gStorage(NL-1, TwoSLdown, ILdown, NL,   TwoSRdown, IRdown);
                              double * LleftBlk  =  Lleft->gStorage(NL-1, TwoSLdown, ILdown, NL,   TwoSLup,   ILup);
                              double * LrightBlk = Lright->gStorage(NL,   TwoSRdown, IRdown, NL+1, TwoSRup,   IRup);
                              
                              char trans = 'T';
                              char notrans = 'N';
                              double alpha = 1.0;
                              double beta = 0.0; //SET
                              
                              dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,LleftBlk,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
               
                              dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,LrightBlk,&dimRdown,&beta,workmem2,&dimLup);
                              
                              int length = dimLup * dimRup;
                              int inc = 1;
                              double value = ddot_(&length, workmem2, &inc, T_up, &inc);

                              const double fact1 = Special::phase( TwoSLup + TwoSRdown + 2 )
                                                 * ( TwoSRup + 1 ) * sqrt( ( TwoSRdown + 1 ) * ( TwoSLup + 1.0 ) )
                                                 * Wigner::wigner6j( TwoSRup, 1, TwoSLup, TwoSLdown, 1, TwoSRdown );
                              const double fact2 = 2 * ( TwoSRup + 1 ) * sqrt( ( TwoSRdown + 1 ) * ( TwoSLup + 1.0 ) )
                                                 * Wigner::wigner6j( TwoSRup, TwoSLdown, 2, 1, 1, TwoSLup )
                                                 * Wigner::wigner6j( TwoSRup, TwoSLdown, 2, 1, 1, TwoSRdown );
                              const int fact3 = (( TwoSRdown == TwoSLup ) ? ( TwoSRup + 1 ) : 0 );

                              d9[0] += fact1  * value;
                              d10[0] += fact2 * value;
                              d11[0] += fact3 * value;

                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
}

double CheMPS2::TwoDM::doD12(TensorT * denT, TensorL * Lleft, TensorL * Lright, double * workmem, double * workmem2, int Irrep_g){

   double d12 = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL, TwoSLup, ILup);
            int dimRup = denBK->gCurrentDim(theindex+1, NL, TwoSLup, ILup);
            if ((dimLup>0) && (dimRup>0)){
            
               int Idown = Irreps::directProd(ILup, Irrep_g);
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                        
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, Idown);
                  int dimRdown = denBK->gCurrentDim(theindex+1, NL+1, TwoSLdown, Idown);
               
                  if ((dimLdown>0) && (dimRdown>0)){
                           
                     double * T_up      =   denT->gStorage(NL,   TwoSLup,   ILup,  NL,   TwoSLup,   ILup);
                     double * T_down    =   denT->gStorage(NL-1, TwoSLdown, Idown, NL+1, TwoSLdown, Idown);
                     double * LleftBlk  =  Lleft->gStorage(NL-1, TwoSLdown, Idown, NL,   TwoSLup,   ILup);
                     double * LrightBlk = Lright->gStorage(NL,   TwoSLup,   ILup,  NL+1, TwoSLdown, Idown);
                              
                     char trans = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta = 0.0; //SET
                              
                     dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,LleftBlk,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
               
                     dgemm_(&notrans,&trans,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,LrightBlk,&dimRup,&beta,workmem2,&dimLup);

                     const double factor = Special::phase( TwoSLdown + 1 - TwoSLup )
                                         * 0.5 * sqrt( ( TwoSLup + 1 ) * ( TwoSLdown + 1.0 ) );

                     int length = dimLup * dimRup;
                     int inc = 1;
                     d12 += factor * ddot_(&length, workmem2, &inc, T_up, &inc);

                  }
               }
            }
         }
      }
   }
   
   return d12;
   
}

double CheMPS2::TwoDM::doD13(TensorT * denT, TensorL * Lleft, TensorS0 * S0right, double * workmem, double * workmem2, int Irrep_g){

   double d13 = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            int dimRup = denBK->gCurrentDim(theindex+1, NL+2, TwoSLup, ILup);
            
            if ((dimLup>0) && (dimRup>0)){
            
               int ILdown = Irreps::directProd(ILup,   Irrep_g);
               int IRdown = Irreps::directProd(ILdown, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                        
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+1, NL,   TwoSLup,   IRdown);
               
                  if ((dimLdown>0) && (dimRdown>0)){
                           
                     double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,   NL+2, TwoSLup, ILup);
                     double * T_down  =    denT->gStorage(NL-1, TwoSLdown, ILdown, NL,   TwoSLup, IRdown);
                     double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, ILdown, NL,   TwoSLup, ILup);
                     double * S0block = S0right->gStorage(NL,   TwoSLup,   IRdown, NL+2, TwoSLup, ILup);
                              
                     char trans = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta = 0.0; //SET
                              
                     dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
               
                     dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,S0block,&dimRdown,&beta,workmem2,&dimLup);
               
                     double factor = -0.5 * sqrt(0.5) * (TwoSLup+1);
                     
                     int length = dimLup * dimRup;
                     int inc = 1;
                     d13 += factor * ddot_(&length, workmem2, &inc, T_up, &inc);

                  }
               }
            }
         }
      }
   }
   
   return d13;
   
}

double CheMPS2::TwoDM::doD14(TensorT * denT, TensorL * Lleft, TensorS0 * S0right, double * workmem, double * workmem2, int Irrep_g){

   double d14 = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            
            if (dimLup>0){
            
               int Idown = Irreps::directProd(ILup, Irrep_g);
               int IRup  = Irreps::directProd(ILup, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  
                  int dimRup   = denBK->gCurrentDim(theindex+1, NL+1, TwoSLdown, IRup);
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, Idown);
                  int dimRdown = denBK->gCurrentDim(theindex+1, NL-1, TwoSLdown, Idown);
               
                  if ((dimLdown>0) && (dimRdown>0) && (dimRup>0)){
                           
                     double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,  NL+1, TwoSLdown, IRup);
                     double * T_down  =    denT->gStorage(NL-1, TwoSLdown, Idown, NL-1, TwoSLdown, Idown);
                     double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, Idown, NL,   TwoSLup,   ILup);
                     double * S0block = S0right->gStorage(NL-1, TwoSLdown, Idown, NL+1, TwoSLdown, IRup);
                              
                     char trans = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta = 0.0; //SET
                              
                     dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
               
                     dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,S0block,&dimRdown,&beta,workmem2,&dimLup);

                     const double factor = Special::phase( TwoSLdown + 1 - TwoSLup ) * 0.5 * sqrt(0.5 * (TwoSLup+1) * (TwoSLdown+1));

                     int length = dimLup * dimRup;
                     int inc = 1;
                     d14 += factor * ddot_(&length, workmem2, &inc, T_up, &inc);

                  }
               }
            }
         }
      }
   }
   
   return d14;
   
}

double CheMPS2::TwoDM::doD15(TensorT * denT, TensorL * Lleft, TensorS1 * S1right, double * workmem, double * workmem2, int Irrep_g){

   double d15 = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            int dimRup = denBK->gCurrentDim(theindex+1, NL+2, TwoSLup, ILup);
            
            if ((dimLup>0) && (dimRup>0)){
            
               int ILdown = Irreps::directProd(ILup,   Irrep_g);
               int IRdown = Irreps::directProd(ILdown, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  for (int TwoSRdown=TwoSLdown-1; TwoSRdown<=TwoSLdown+1; TwoSRdown+=2){
                        
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+1, NL,   TwoSRdown, IRdown);
               
                     if ((dimLdown>0) && (dimRdown>0)){
                           
                        double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,   NL+2, TwoSLup,   ILup);
                        double * T_down  =    denT->gStorage(NL-1, TwoSLdown, ILdown, NL,   TwoSRdown, IRdown);
                        double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, ILdown, NL,   TwoSLup,   ILup);
                        double * S1block = S1right->gStorage(NL,   TwoSRdown, IRdown, NL+2, TwoSLup,   ILup);
                              
                        char trans = 'T';
                        char notrans = 'N';
                        double alpha = 1.0;
                        double beta = 0.0; //SET
                              
                        dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
               
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,S1block,&dimRdown,&beta,workmem2,&dimLup);
               
                        const double factor = Special::phase( TwoSLdown + TwoSLup + 1 )
                                            * ( TwoSLup + 1 ) * sqrt((TwoSRdown+1)/3.0) * Wigner::wigner6j( 1, 1, 2, TwoSLup, TwoSRdown, TwoSLdown );
                        
                        int length = dimLup * dimRup;
                        int inc = 1;
                        d15 += factor * ddot_(&length, workmem2, &inc, T_up, &inc);
                     
                     }
                  }
               }
            }
         }
      }
   }
   
   return d15;
   
}

double CheMPS2::TwoDM::doD16(TensorT * denT, TensorL * Lleft, TensorS1 * S1right, double * workmem, double * workmem2, int Irrep_g){

   double d16 = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            
            if (dimLup>0){
            
               int Idown = Irreps::directProd(ILup, Irrep_g);
               int IRup  = Irreps::directProd(ILup, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  for (int TwoSRup=TwoSLup-1; TwoSRup<=TwoSLup+1; TwoSRup+=2){
                  
                     int dimRup   = denBK->gCurrentDim(theindex+1, NL+1, TwoSRup,   IRup);
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, Idown);
                     int dimRdown = denBK->gCurrentDim(theindex+1, NL-1, TwoSLdown, Idown);
                  
                     if ((dimLdown>0) && (dimRdown>0) && (dimRup>0)){
                              
                        double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,  NL+1, TwoSRup,   IRup);
                        double * T_down  =    denT->gStorage(NL-1, TwoSLdown, Idown, NL-1, TwoSLdown, Idown);
                        double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, Idown, NL,   TwoSLup,   ILup);
                        double * S1block = S1right->gStorage(NL-1, TwoSLdown, Idown, NL+1, TwoSRup,   IRup);
                                 
                        char trans = 'T';
                        char notrans = 'N';
                        double alpha = 1.0;
                        double beta = 0.0; //SET
                                 
                        dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
                  
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,S1block,&dimRdown,&beta,workmem2,&dimLup);

                        const double factor = Special::phase( TwoSRup + TwoSLdown + 2 ) * (TwoSRup+1) * sqrt((TwoSLup+1)/3.0)
                                            * Wigner::wigner6j( 1, 1, 2, TwoSRup, TwoSLdown, TwoSLup );

                        int length = dimLup * dimRup;
                        int inc = 1;
                        d16 += factor * ddot_(&length, workmem2, &inc, T_up, &inc);
                     
                     }
                  }
               }
            }
         }
      }
   }
   
   return d16;
   
}

double CheMPS2::TwoDM::doD17orD21(TensorT * denT, TensorL * Lleft, TensorF0 * F0right, double * workmem, double * workmem2, int Irrep_g, bool shouldIdoD17){

   double total = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            
            if (dimLup>0){
            
               int ILdown = Irreps::directProd(ILup,   Irrep_g);
               int IRdown = Irreps::directProd(ILdown, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  
                  int dimRup   = denBK->gCurrentDim(theindex+1, NL,   TwoSLup,   ILup);
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+1, NL,   TwoSLup,   IRdown);
               
                  if ((dimLdown>0) && (dimRdown>0) && (dimRup>0)){
                           
                     double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,   NL, TwoSLup, ILup);
                     double * T_down  =    denT->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup, IRdown);
                     double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup, ILup);
                     double * F0block = (shouldIdoD17) ? F0right->gStorage(NL,   TwoSLup,   IRdown, NL, TwoSLup, ILup)
                                                       : F0right->gStorage(NL,   TwoSLup,   ILup,   NL, TwoSLup, IRdown) ;
                              
                     char trans = 'T';
                     char notrans = 'N';
                     char var = (shouldIdoD17) ? notrans : trans;
                     double alpha = 1.0;
                     double beta = 0.0; //SET
                     int dimvar = (shouldIdoD17) ? dimRdown : dimRup ;
                              
                     dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
               
                     dgemm_(&notrans,&var,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,F0block,&dimvar,&beta,workmem2,&dimLup);
                     
                     int length = dimLup * dimRup;
                     int inc = 1;
                     total += sqrt(0.5) * 0.5 * (TwoSLup+1) * ddot_(&length, workmem2, &inc, T_up, &inc);

                  }
               }
            }
         }
      }
   }
   
   return total;
   
}

double CheMPS2::TwoDM::doD18orD22(TensorT * denT, TensorL * Lleft, TensorF0 * F0right, double * workmem, double * workmem2, int Irrep_g, bool shouldIdoD18){

   double total = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            
            if (dimLup>0){
            
               int Idown = Irreps::directProd(ILup, Irrep_g);
               int IRup  = Irreps::directProd(ILup, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  
                  int dimRup   = denBK->gCurrentDim(theindex+1, NL+1, TwoSLdown, IRup);
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, Idown);
                  int dimRdown = denBK->gCurrentDim(theindex+1, NL+1, TwoSLdown, Idown);
               
                  if ((dimLdown>0) && (dimRdown>0) && (dimRup>0)){
                           
                     double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,  NL+1, TwoSLdown, IRup);
                     double * T_down  =    denT->gStorage(NL-1, TwoSLdown, Idown, NL+1, TwoSLdown, Idown);
                     double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, Idown, NL, TwoSLup, ILup);
                     double * F0block = (shouldIdoD18) ? F0right->gStorage(NL+1, TwoSLdown, Idown, NL+1, TwoSLdown, IRup)
                                                       : F0right->gStorage(NL+1, TwoSLdown, IRup,  NL+1, TwoSLdown, Idown) ;
                              
                     char trans = 'T';
                     char notrans = 'N';
                     char var = (shouldIdoD18) ? notrans : trans;
                     double alpha = 1.0;
                     double beta = 0.0; //SET
                     int dimvar = (shouldIdoD18) ? dimRdown : dimRup ;
                              
                     dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
               
                     dgemm_(&notrans,&var,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,F0block,&dimvar,&beta,workmem2,&dimLup);

                     const double factor = Special::phase( TwoSLdown + 1 - TwoSLup ) * 0.5 * sqrt(0.5*(TwoSLup+1)*(TwoSLdown+1));

                     int length = dimLup * dimRup;
                     int inc = 1;
                     total += factor * ddot_(&length, workmem2, &inc, T_up, &inc);

                  }
               }
            }
         }
      }
   }
   
   return total;
   
}

double CheMPS2::TwoDM::doD19orD23(TensorT * denT, TensorL * Lleft, TensorF1 * F1right, double * workmem, double * workmem2, int Irrep_g, bool shouldIdoD19){

   double total = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex, NL, TwoSLup, ILup);
            
            if (dimLup>0){
            
               int ILdown = Irreps::directProd(ILup,   Irrep_g);
               int IRdown = Irreps::directProd(ILdown, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  for (int TwoSRdown=TwoSLdown-1; TwoSRdown<=TwoSLdown+1; TwoSRdown+=2){
                  
                     int dimRup   = denBK->gCurrentDim(theindex+1, NL,   TwoSLup,   ILup);
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+1, NL,   TwoSRdown, IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0) && (dimRup>0)){
                              
                        double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,   NL, TwoSLup,   ILup);
                        double * T_down  =    denT->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSRdown, IRdown);
                        double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup,   ILup);
                        double * F1block = (shouldIdoD19) ? F1right->gStorage(NL, TwoSRdown, IRdown, NL, TwoSLup,   ILup)
                                                          : F1right->gStorage(NL, TwoSLup,   ILup,   NL, TwoSRdown, IRdown) ;
                                 
                        char trans = 'T';
                        char notrans = 'N';
                        char var = (shouldIdoD19) ? notrans : trans;
                        double alpha = 1.0;
                        double beta = 0.0; //SET
                        int dimvar = (shouldIdoD19) ? dimRdown : dimRup ;
                                 
                        dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
                  
                        dgemm_(&notrans,&var,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,F1block,&dimvar,&beta,workmem2,&dimLup);
                  
                        double factor = 0.0;
                        if (shouldIdoD19){
                           factor = Special::phase( TwoSLdown + TwoSRdown - 1 )
                                  * ( TwoSRdown + 1 ) * sqrt((TwoSLup+1)/3.0) * Wigner::wigner6j( 1, 1, 2, TwoSLup, TwoSRdown, TwoSLdown );
                        } else {
                           factor = Special::phase( TwoSLdown + TwoSLup - 1 )
                                  * ( TwoSLup + 1 ) * sqrt((TwoSRdown+1)/3.0) * Wigner::wigner6j( 1, 1, 2, TwoSLup, TwoSRdown, TwoSLdown );
                        }
                        
                        int length = dimLup * dimRup;
                        int inc = 1;
                        total += factor * ddot_(&length,workmem2,&inc,T_up,&inc);
                     
                     }
                  }
               }
            }
         }
      }
   }
   
   return total;
   
}

double CheMPS2::TwoDM::doD20orD24(TensorT * denT, TensorL * Lleft, TensorF1 * F1right, double * workmem, double * workmem2, int Irrep_g, bool shouldIdoD20){

   double total = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            
            if (dimLup>0){
            
               int Idown = Irreps::directProd(ILup, Irrep_g);
               int IRup  = Irreps::directProd(ILup, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  for (int TwoSRup=TwoSLup-1; TwoSRup<=TwoSLup+1; TwoSRup+=2){
                  
                     int dimRup   = denBK->gCurrentDim(theindex+1, NL+1, TwoSRup,   IRup);
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, Idown);
                     int dimRdown = denBK->gCurrentDim(theindex+1, NL+1, TwoSLdown, Idown);
                  
                     if ((dimLdown>0) && (dimRdown>0) && (dimRup>0)){
                              
                        double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,  NL+1, TwoSRup,   IRup);
                        double * T_down  =    denT->gStorage(NL-1, TwoSLdown, Idown, NL+1, TwoSLdown, Idown);
                        double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, Idown, NL,   TwoSLup,   ILup);
                        double * F1block = (shouldIdoD20) ? F1right->gStorage(NL+1, TwoSLdown, Idown, NL+1, TwoSRup,   IRup)
                                                          : F1right->gStorage(NL+1, TwoSRup,   IRup,  NL+1, TwoSLdown, Idown) ;
                                 
                        char trans = 'T';
                        char notrans = 'N';
                        char var = (shouldIdoD20) ? notrans : trans;
                        double alpha = 1.0;
                        double beta = 0.0; //SET
                        int dimvar = (shouldIdoD20) ? dimRdown : dimRup ;
                                 
                        dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
                  
                        dgemm_(&notrans,&var,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,F1block,&dimvar,&beta,workmem2,&dimLup);
                  
                        double factor = 0.0;
                        if (shouldIdoD20){
                           factor = Special::phase( 2 * TwoSLup )
                                  * sqrt((TwoSLup+1)*(TwoSRup+1)*(TwoSLdown+1)/3.0) * Wigner::wigner6j( 1, 1, 2, TwoSRup, TwoSLdown, TwoSLup );
                        } else {
                           factor = Special::phase( 2 * TwoSLup + TwoSRup - TwoSLdown )
                                  * (TwoSRup+1) * sqrt((TwoSLup+1)/3.0) * Wigner::wigner6j( 1, 1, 2, TwoSRup, TwoSLdown, TwoSLup );
                        }
                        
                        int length = dimLup * dimRup;
                        int inc = 1;
                        total += factor * ddot_(&length, workmem2, &inc, T_up, &inc);
                     
                     }
                  }
               }
            }
         }
      }
   }
   
   return total;
   
}




