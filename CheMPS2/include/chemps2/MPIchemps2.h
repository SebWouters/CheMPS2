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

#ifndef MPICHEMPS2_CHEMPS2_H
#define MPICHEMPS2_CHEMPS2_H

//#define CHEMPS2_MPI_COMPILATION

#ifdef CHEMPS2_MPI_COMPILATION

   #include <mpi.h>
   #include <assert.h>
   #include "Tensor.h"

   #define MPI_CHEMPS2_MASTER   0

   //Assign diagrams which are not specifically owned round-robin

   #define MPI_CHEMPS2_4D1AB    1
   #define MPI_CHEMPS2_4D2AB    2
   #define MPI_CHEMPS2_4I1AB    3
   #define MPI_CHEMPS2_4I2AB    4
   #define MPI_CHEMPS2_4F1AB    5
   #define MPI_CHEMPS2_4F2AB    6
   #define MPI_CHEMPS2_4G1AB    7
   #define MPI_CHEMPS2_4G2AB    8

   #define MPI_CHEMPS2_4D3ABCD  9
   #define MPI_CHEMPS2_4D4ABCD  10
   #define MPI_CHEMPS2_4I3ABCD  11
   #define MPI_CHEMPS2_4I4ABCD  12
   #define MPI_CHEMPS2_4F3ABCD  13
   #define MPI_CHEMPS2_4F4ABCD  14
   #define MPI_CHEMPS2_4G3ABCD  15
   #define MPI_CHEMPS2_4G4ABCD  16

   #define MPI_CHEMPS2_4E1      17
   #define MPI_CHEMPS2_4E2      18
   #define MPI_CHEMPS2_4H1      19
   #define MPI_CHEMPS2_4H2      20

   #define MPI_CHEMPS2_4E3A     21
   #define MPI_CHEMPS2_4E3B     22
   #define MPI_CHEMPS2_4E4A     23
   #define MPI_CHEMPS2_4E4B     24
   #define MPI_CHEMPS2_4H3A     25
   #define MPI_CHEMPS2_4H3B     26
   #define MPI_CHEMPS2_4H4A     27
   #define MPI_CHEMPS2_4H4B     28

   #define MPI_CHEMPS2_5A1      29
   #define MPI_CHEMPS2_5A2      30
   #define MPI_CHEMPS2_5A3      31
   #define MPI_CHEMPS2_5A4      32

   #define MPI_CHEMPS2_5B1      33
   #define MPI_CHEMPS2_5B2      34
   #define MPI_CHEMPS2_5B3      35
   #define MPI_CHEMPS2_5B4      36

   #define MPI_CHEMPS2_5C1      37
   #define MPI_CHEMPS2_5C2      38
   #define MPI_CHEMPS2_5C3      39
   #define MPI_CHEMPS2_5C4      40

   #define MPI_CHEMPS2_5D1      41
   #define MPI_CHEMPS2_5D2      42
   #define MPI_CHEMPS2_5D3      43
   #define MPI_CHEMPS2_5D4      44

   #define MPI_CHEMPS2_5E1      45
   #define MPI_CHEMPS2_5E2      46
   #define MPI_CHEMPS2_5E3      47
   #define MPI_CHEMPS2_5E4      48

   #define MPI_CHEMPS2_5F1      49
   #define MPI_CHEMPS2_5F2      50
   #define MPI_CHEMPS2_5F3      51
   #define MPI_CHEMPS2_5F4      52

   #define MPI_CHEMPS2_OFFSET   53

#endif

namespace CheMPS2{
/** MPIchemps2 class
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date May 27, 2015
    
    The MPIchemps2 class contains the bookkeeping for MPI */
   class MPIchemps2{

      public:
      
         //! Constructor
         MPIchemps2(){}
         
         //! Destructor
         virtual ~MPIchemps2(){}
         
         //! Get the number of MPI processes
         /** \return The number of MPI processes */
         static int mpi_size(){
            #ifdef CHEMPS2_MPI_COMPILATION
               int size;
               MPI_Comm_size( MPI_COMM_WORLD, &size );
               return size;
            #else
               return 1;
            #endif
         }

         //! Get the rank of this MPI process
         /** \return The rank of this MPI process */
         static int mpi_rank(){
            #ifdef CHEMPS2_MPI_COMPILATION
               int rank;
               MPI_Comm_rank( MPI_COMM_WORLD, &rank );
               return rank;
            #else
               return 0;
            #endif
         }
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Initialize MPI
         static void mpi_init(){
            int zero = 0;
            MPI_Init( &zero, NULL );
         }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Finalize MPI
         static void mpi_finalize(){
            MPI_Finalize();
         }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Get the owner of the X-tensors
         static int owner_x(){ return MPI_CHEMPS2_MASTER; }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Get the owner of a certain {A,B,Sigma0,Sigma1}-tensor
         /** \param index1 The first  DMRG lattice index of the tensor
             \param index2 The second DMRG lattice index of the tensor
             \return The owner rank */
         static int owner_absigma(const int index1, const int index2){ // 1 <= proc < 1 + L*(L+1)/2
            assert( index1 <= index2 );
            return ( 1 + index1 + (index2*(index2+1))/2 ) % mpi_size();
         }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Get the owner of a certain {C,D,F0,F1}-tensor
         /** \param L The number of active space orbitals
             \param index1 The first  DMRG lattice index of the tensor
             \param index2 The second DMRG lattice index of the tensor
             \return The owner rank */
         static int owner_cdf(const int L, const int index1, const int index2){ // 1 + L*(L+1)/2 <= proc < 1 + L*(L+1)
            assert( index1 <= index2 );
            return ( 1 + (L*(L+1))/2 + index1 + (index2*(index2+1))/2 ) % mpi_size();
         }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Get the owner of a certain Q-tensor
         /** \param L The number of active space orbitals
             \param index The DMRG lattice index of the tensor
             \return The owner rank */
         static int owner_q(const int L, const int index){ // 1 + L*(L+1) <= proc < 1 + L*(L+2)
            return ( 1 + L*(L+1) + index ) % mpi_size();
         }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Get the owner of a certain 3-index tensor for the 3-RDM
         /** \param L The number of active space orbitals
             \param index1 The first  DMRG lattice index of the tensor
             \param index2 The second DMRG lattice index of the tensor
             \param index3 The third  DMRG lattice index of the tensor
             \return The owner rank */
         static int owner_3rdm_diagram(const int L, const int index1, const int index2, const int index3){ // 1 + L*(L+1) <= proc < 1 + L*(L+1) + L*(L+1)*(L+2)/6
            assert( index1 <= index2 );
            assert( index2 <= index3 );
            return ( 1 + L*(L+1) + index1 + (index2*(index2+1))/2 + (index3*(index3+1)*(index3+2))/6 ) % mpi_size();
         }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Get the owner of the 1c, 1d, 2d, 3e, and 3h diagrams of the effective Hamiltonian
         /** \return The owner rank */
         static int owner_1cd2d3eh(){ return MPI_CHEMPS2_MASTER; }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Get the owner of a specific diagram (or group of diagrams) of the effective Hamiltonian
         /** \param L The number of active space orbitals
             \param macro The macro number of the specific diagram (or group of diagrams) as defined in the file MPIchemps2.h
             \return The owner rank */
         static int owner_specific_diagram(const int L, const int macro){
            return (macro + L*(L+2)) % mpi_size();
         }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Get the owner of a specific excitation
         /** \param L The number of active space orbitals
             \param excitation The number of the specific excitation
             \return The owner rank */
         static int owner_specific_excitation(const int L, const int excitation){
            return (MPI_CHEMPS2_OFFSET + L*(L+2) + excitation) % mpi_size();
         }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Broadcast a tensor
         /** \param object The tensor to be broadcasted
             \param ROOT The MPI process which should broadcast */
         static void broadcast_tensor(Tensor * object, int ROOT){
            int arraysize = object->gKappa2index(object->gNKappa());
            MPI_Bcast(object->gStorage(), arraysize, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
         }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Broadcast an array of doubles
         /** \param array The array to be broadcasted
             \param length The length of the array
             \param ROOT The MPI process which should broadcast */
         static void broadcast_array_double(double * array, int length, int ROOT){
            MPI_Bcast(array, length, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
         }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Broadcast an array of integers
         /** \param array The array to be broadcasted
             \param length The length of the array
             \param ROOT The MPI process which should broadcast */
         static void broadcast_array_int(int * array, int length, int ROOT){
            MPI_Bcast(array, length, MPI_INT, ROOT, MPI_COMM_WORLD);
         }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Check whether all processes agree on a boolean
         /** \param mybool The process's boolean
             \return Whether the booleans of all processes are equal */
         static bool all_booleans_equal(const bool mybool){
            int my_value = ( mybool ) ? 1 : 0 ;
            int tot_value;
            MPI_Allreduce(&my_value, &tot_value, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            return ( my_value * MPIchemps2::mpi_size() == tot_value ); // Only true if mybool is the same for all processes
         }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Send a tensor from one process to another
         /** \param object The tensor to be sent
             \param SENDER The MPI process which should send the tensor
             \param RECEIVER The MPI process which should receive the tensor
             \param tag A tag which should be the same for the sender and receiver to make sure that the communication is desired */
         static void sendreceive_tensor(Tensor * object, int SENDER, int RECEIVER, int tag){
            if ( SENDER != RECEIVER ){
               const int MPIRANK = mpi_rank();
               if ( SENDER == MPIRANK ){
                  int arraysize = object->gKappa2index(object->gNKappa());
                  MPI_Send(object->gStorage(), arraysize, MPI_DOUBLE, RECEIVER, tag, MPI_COMM_WORLD);
               }
               if ( RECEIVER == MPIRANK ){
                  int arraysize = object->gKappa2index(object->gNKappa());
                  MPI_Recv(object->gStorage(), arraysize, MPI_DOUBLE, SENDER, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               }
            }
         }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Add arrays of all processes and give result to ROOT
         /** \param vec_in The array which should be added
             \param vec_out The array where the result should be stored
             \param size The size of the array
             \param ROOT The MPI process which should have the result vector */
         static void reduce_array_double(double * vec_in, double * vec_out, int size, int ROOT){
            MPI_Reduce(vec_in, vec_out, size, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
         }
         #endif
         
         #ifdef CHEMPS2_MPI_COMPILATION
         //! Add arrays of all processes and give everyone the result
         /** \param vec_in The array which should be added
             \param vec_out The array where the result should be stored
             \param size The size of the array */
         static void allreduce_array_double(double * vec_in, double * vec_out, int size){
            MPI_Allreduce(vec_in, vec_out, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         }
         #endif

   };
}

#endif

