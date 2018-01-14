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

#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>

#include "Sobject.h"
#include "TensorT.h"
#include "SyBookkeeper.h"
#include "Lapack.h"
#include "MPIchemps2.h"
#include "Wigner.h"
#include "Special.h"

using std::min;
using std::max;

CheMPS2::Sobject::Sobject( const int index, SyBookkeeper * denBK ){

   this->index = index;
   this->denBK = denBK;
   Ilocal1 = denBK->gIrrep( index     );
   Ilocal2 = denBK->gIrrep( index + 1 );

   nKappa = 0;

   for ( int NL = denBK->gNmin( index ); NL <= denBK->gNmax( index ); NL++ ){
      for ( int TwoSL = denBK->gTwoSmin( index, NL ); TwoSL <= denBK->gTwoSmax( index, NL ); TwoSL += 2 ){
         for ( int IL = 0; IL < denBK->getNumberOfIrreps(); IL++ ){
            const int dimL = denBK->gCurrentDim( index, NL, TwoSL, IL );
            if ( dimL > 0 ){
               for ( int N1 = 0; N1 <= 2; N1++ ){
                  for ( int N2 = 0; N2 <= 2; N2++ ){
                     const int NR = NL + N1 + N2;
                     const int IM = (( N1 == 1 ) ? Irreps::directProd( IL, Ilocal1 ) : IL );
                     const int IR = (( N2 == 1 ) ? Irreps::directProd( IM, Ilocal2 ) : IM );
                     const int TwoJmin = ( N1 + N2 ) % 2;
                     const int TwoJmax = ((( N1 == 1 ) && ( N2 == 1 )) ? 2 : TwoJmin );
                     for ( int TwoJ = TwoJmin; TwoJ <= TwoJmax; TwoJ += 2 ){
                        for ( int TwoSR = TwoSL - TwoJ; TwoSR <= TwoSL + TwoJ; TwoSR += 2 ){
                           if ( TwoSR >= 0 ){
                              const int dimR = denBK->gCurrentDim( index + 2, NR, TwoSR, IR );
                              if ( dimR > 0 ){
                                 nKappa++;
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

   sectorNL    = new int[ nKappa ];
   sectorTwoSL = new int[ nKappa ];
   sectorIL    = new int[ nKappa ];
   sectorN1    = new int[ nKappa ];
   sectorN2    = new int[ nKappa ];
   sectorTwoJ  = new int[ nKappa ];
   sectorNR    = new int[ nKappa ];
   sectorTwoSR = new int[ nKappa ];
   sectorIR    = new int[ nKappa ];
   kappa2index = new int[ nKappa + 1 ];
   kappa2index[ 0 ] = 0;

   nKappa = 0;

   for ( int NL = denBK->gNmin( index ); NL <= denBK->gNmax( index ); NL++ ){
      for ( int TwoSL = denBK->gTwoSmin( index, NL ); TwoSL <= denBK->gTwoSmax( index, NL ); TwoSL += 2 ){
         for ( int IL = 0; IL < denBK->getNumberOfIrreps(); IL++ ){
            const int dimL = denBK->gCurrentDim( index, NL, TwoSL, IL );
            if ( dimL > 0 ){
               for ( int N1 = 0; N1 <= 2; N1++ ){
                  for ( int N2 = 0; N2 <= 2; N2++ ){
                     const int NR = NL + N1 + N2;
                     const int IM = (( N1 == 1 ) ? Irreps::directProd( IL, Ilocal1 ) : IL );
                     const int IR = (( N2 == 1 ) ? Irreps::directProd( IM, Ilocal2 ) : IM );
                     const int TwoJmin = ( N1 + N2 ) % 2;
                     const int TwoJmax = ((( N1 == 1 ) && ( N2 == 1 )) ? 2 : TwoJmin );
                     for ( int TwoJ = TwoJmin; TwoJ <= TwoJmax; TwoJ += 2 ){
                        for ( int TwoSR = TwoSL - TwoJ; TwoSR <= TwoSL + TwoJ; TwoSR += 2 ){
                           if ( TwoSR >= 0 ){
                              const int dimR = denBK->gCurrentDim( index + 2, NR, TwoSR, IR );
                              if ( dimR > 0 ){
                                 sectorNL   [ nKappa ] = NL;
                                 sectorTwoSL[ nKappa ] = TwoSL;
                                 sectorIL   [ nKappa ] = IL;
                                 sectorN1   [ nKappa ] = N1;
                                 sectorN2   [ nKappa ] = N2;
                                 sectorTwoJ [ nKappa ] = TwoJ;
                                 sectorNR   [ nKappa ] = NR;
                                 sectorTwoSR[ nKappa ] = TwoSR;
                                 sectorIR   [ nKappa ] = IR;
                                 nKappa++;
                                 kappa2index[ nKappa ] = kappa2index[ nKappa - 1 ] + dimL * dimR;
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

   storage = new double[ kappa2index[ nKappa ] ];

   reorder = new int[ nKappa ];
   for ( int cnt = 0; cnt < nKappa; cnt++ ){ reorder[ cnt ] = cnt; }
   bool sorted = false;
   while ( sorted == false ){ //Bubble sort so that blocksize(reorder[i]) >= blocksize(reorder[i+1]), with blocksize(k) = kappa2index[k+1]-kappa2index[k]
      sorted = true;
      for ( int cnt = 0; cnt < nKappa - 1; cnt++ ){
         const int index1 = reorder[ cnt ];
         const int index2 = reorder[ cnt + 1 ];
         const int size1  = kappa2index[ index1 + 1 ] - kappa2index[ index1 ];
         const int size2  = kappa2index[ index2 + 1 ] - kappa2index[ index2 ];
         if ( size1 < size2 ){
            sorted = false;
            reorder[ cnt ]     = index2;
            reorder[ cnt + 1 ] = index1;
         }
      }
   }

}

CheMPS2::Sobject::~Sobject(){

   delete [] sectorNL;
   delete [] sectorTwoSL;
   delete [] sectorIL;
   delete [] sectorN1;
   delete [] sectorN2;
   delete [] sectorTwoJ;
   delete [] sectorNR;
   delete [] sectorTwoSR;
   delete [] sectorIR;
   delete [] kappa2index;
   delete [] storage;
   delete [] reorder;

}

int CheMPS2::Sobject::gNKappa() const { return nKappa; }

double * CheMPS2::Sobject::gStorage() { return storage; }

int CheMPS2::Sobject::gReorder( const int ikappa ) const{ return reorder[ ikappa ]; }

int CheMPS2::Sobject::gKappa( const int NL, const int TwoSL, const int IL, const int N1, const int N2, const int TwoJ, const int NR, const int TwoSR, const int IR ) const{

   for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){
      if (( sectorNL   [ ikappa ] == NL    ) &&
          ( sectorTwoSL[ ikappa ] == TwoSL ) &&
          ( sectorIL   [ ikappa ] == IL    ) &&
          ( sectorN1   [ ikappa ] == N1    ) &&
          ( sectorN2   [ ikappa ] == N2    ) &&
          ( sectorTwoJ [ ikappa ] == TwoJ  ) &&
          ( sectorNR   [ ikappa ] == NR    ) &&
          ( sectorTwoSR[ ikappa ] == TwoSR ) &&
          ( sectorIR   [ ikappa ] == IR    )){ return ikappa; }
   }

   return -1;

}

int CheMPS2::Sobject::gKappa2index( const int kappa ) const{ return kappa2index[ kappa ]; }

double * CheMPS2::Sobject::gStorage( const int NL, const int TwoSL, const int IL, const int N1, const int N2, const int TwoJ, const int NR, const int TwoSR, const int IR ){

   const int kappa = gKappa( NL, TwoSL, IL, N1, N2, TwoJ, NR, TwoSR, IR );
   if ( kappa == -1 ){ return NULL; }
   return storage + kappa2index[ kappa ];

}

int CheMPS2::Sobject::gIndex() const { return index; }

int CheMPS2::Sobject::gNL   ( const int ikappa ) const{ return sectorNL   [ ikappa ]; }
int CheMPS2::Sobject::gTwoSL( const int ikappa ) const{ return sectorTwoSL[ ikappa ]; }
int CheMPS2::Sobject::gIL   ( const int ikappa ) const{ return sectorIL   [ ikappa ]; }
int CheMPS2::Sobject::gN1   ( const int ikappa ) const{ return sectorN1   [ ikappa ]; }
int CheMPS2::Sobject::gN2   ( const int ikappa ) const{ return sectorN2   [ ikappa ]; }
int CheMPS2::Sobject::gTwoJ ( const int ikappa ) const{ return sectorTwoJ [ ikappa ]; }
int CheMPS2::Sobject::gNR   ( const int ikappa ) const{ return sectorNR   [ ikappa ]; }
int CheMPS2::Sobject::gTwoSR( const int ikappa ) const{ return sectorTwoSR[ ikappa ]; }
int CheMPS2::Sobject::gIR   ( const int ikappa ) const{ return sectorIR   [ ikappa ]; }

void CheMPS2::Sobject::Join( TensorT * Tleft, TensorT * Tright ){

   #pragma omp parallel for schedule(dynamic)
   for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){

      const int NL    = sectorNL   [ ikappa ];
      const int TwoSL = sectorTwoSL[ ikappa ];
      const int IL    = sectorIL   [ ikappa ];

      const int NR    = sectorNR   [ ikappa ];
      const int TwoSR = sectorTwoSR[ ikappa ];
      const int IR    = sectorIR   [ ikappa ];

      const int TwoJ  = sectorTwoJ [ ikappa ];
      const int N1    = sectorN1   [ ikappa ];
      const int N2    = sectorN2   [ ikappa ];
      const int TwoS1 = (( N1 == 1 ) ? 1 : 0 );
      const int TwoS2 = (( N2 == 1 ) ? 1 : 0 );
      const int fase  = Special::phase( TwoSL + TwoSR + TwoS1 + TwoS2 );

      // Clear block
      int dimL = denBK->gCurrentDim( index,     NL, TwoSL, IL ); // dimL > 0, checked at creation
      int dimR = denBK->gCurrentDim( index + 2, NR, TwoSR, IR ); // dimR > 0, checked at creation
      double * block_s = storage + kappa2index[ ikappa ];
      for ( int cnt = 0; cnt < dimL * dimR; cnt++ ){ block_s[ cnt ] = 0.0; }

      // Central symmetry sectors
      const int NM = NL + N1;
      const int IM = (( TwoS1 == 1 ) ? Irreps::directProd( IL, Ilocal1 ) : IL );
      const int TwoJMlower = max( abs( TwoSL - TwoS1 ), abs( TwoSR - TwoS2 ) );
      const int TwoJMupper = min(    ( TwoSL + TwoS1 ),    ( TwoSR + TwoS2 ) );
      for ( int TwoJM = TwoJMlower; TwoJM <= TwoJMupper; TwoJM += 2 ){
         int dimM = denBK->gCurrentDim( index + 1, NM, TwoJM, IM );
         if ( dimM > 0 ){
            double * block_left  = Tleft ->gStorage( NL, TwoSL, IL, NM, TwoJM, IM );
            double * block_right = Tright->gStorage( NM, TwoJM, IM, NR, TwoSR, IR );
            double prefactor = fase
                             * sqrt( 1.0 * ( TwoJ + 1 ) * ( TwoJM + 1 ) )
                             * Wigner::wigner6j( TwoSL, TwoSR, TwoJ, TwoS2, TwoS1, TwoJM );
            char notrans = 'N';
            double add = 1.0;
            dgemm_( &notrans, &notrans, &dimL, &dimR, &dimM, &prefactor, block_left, &dimL, block_right, &dimM, &add, block_s, &dimL );
         }
      }
   }

}

double CheMPS2::Sobject::Split( TensorT * Tleft, TensorT * Tright, const int virtualdimensionD, const bool movingright, const bool change ){

   #ifdef CHEMPS2_MPI_COMPILATION
   const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #endif

   // Get the number of central sectors
   int nCenterSectors = 0;
   for ( int NM = denBK->gNmin( index + 1 ); NM <= denBK->gNmax( index + 1 ); NM++ ){
      for ( int TwoSM = denBK->gTwoSmin( index + 1, NM ); TwoSM <= denBK->gTwoSmax( index + 1, NM ); TwoSM += 2 ){
         for ( int IM = 0; IM < denBK->getNumberOfIrreps(); IM++ ){
            const int dimM = denBK->gFCIdim( index + 1, NM, TwoSM, IM ); //FCIdim !! Whether possible hence.
            if ( dimM > 0 ){
               nCenterSectors++;
            }
         }
      }
   }

   // Get the labels of the central sectors
   int * SplitSectNM    = new int[ nCenterSectors ];
   int * SplitSectTwoJM = new int[ nCenterSectors ];
   int * SplitSectIM    = new int[ nCenterSectors ];
   nCenterSectors = 0;
   for ( int NM = denBK->gNmin( index + 1 ); NM <= denBK->gNmax( index + 1 ); NM++ ){
      for ( int TwoSM = denBK->gTwoSmin( index + 1, NM ); TwoSM <= denBK->gTwoSmax( index + 1, NM ); TwoSM += 2 ){
         for ( int IM = 0; IM < denBK->getNumberOfIrreps(); IM++ ){
            const int dimM = denBK->gFCIdim( index + 1, NM, TwoSM, IM ); //FCIdim !! Whether possible hence.
            if ( dimM > 0 ){
               SplitSectNM   [ nCenterSectors ] = NM;
               SplitSectTwoJM[ nCenterSectors ] = TwoSM;
               SplitSectIM   [ nCenterSectors ] = IM;
               nCenterSectors++;
            }
         }
      }
   }

   // Only MPI_CHEMPS2_MASTER performs SVD --> Allocate memory
   double ** Lambdas = NULL;
   double ** Us      = NULL;
   double ** VTs     = NULL;
   int * CenterDims  = NULL;
   int * DimLtotal   = NULL;
   int * DimRtotal   = NULL;

   #ifdef CHEMPS2_MPI_COMPILATION
   if ( am_i_master ){
   #endif

   Lambdas = new double*[ nCenterSectors ];
   Us      = new double*[ nCenterSectors ];
   VTs     = new double*[ nCenterSectors ];
   CenterDims  = new int[ nCenterSectors ];
   DimLtotal   = new int[ nCenterSectors ];
   DimRtotal   = new int[ nCenterSectors ];

   //PARALLEL
   #pragma omp parallel for schedule(dynamic)
   for ( int iCenter = 0; iCenter < nCenterSectors; iCenter++ ){

      //Determine left and right dimensions contributing to the center block iCenter
      DimLtotal[ iCenter ] = 0;
      for ( int NL = SplitSectNM[ iCenter ] - 2; NL <= SplitSectNM[ iCenter ]; NL++ ){
         const int TwoS1 = (( NL + 1 == SplitSectNM[ iCenter ] ) ? 1 : 0 );
         for ( int TwoSL = SplitSectTwoJM[ iCenter ] - TwoS1; TwoSL <= SplitSectTwoJM[ iCenter ] + TwoS1; TwoSL += 2 ){
            if ( TwoSL >= 0 ){
               const int IL = (( TwoS1 == 1 ) ? Irreps::directProd( Ilocal1, SplitSectIM[ iCenter ] ) : SplitSectIM[ iCenter ] );
               const int dimL = denBK->gCurrentDim( index, NL, TwoSL, IL );
               if ( dimL > 0 ){
                  DimLtotal[ iCenter ] += dimL;
               }
            }
         }
      }
      DimRtotal[ iCenter ] = 0;
      for ( int NR = SplitSectNM[ iCenter ]; NR <= SplitSectNM[ iCenter ] + 2; NR++ ){
         const int TwoS2 = (( NR == SplitSectNM[ iCenter ] + 1 ) ? 1 : 0 );
         for ( int TwoSR = SplitSectTwoJM[ iCenter ] - TwoS2; TwoSR <= SplitSectTwoJM[ iCenter ] + TwoS2; TwoSR += 2 ){
            if ( TwoSR >= 0 ){
               const int IR = (( TwoS2 == 1 ) ? Irreps::directProd( Ilocal2, SplitSectIM[ iCenter ] ) : SplitSectIM[ iCenter ] );
               const int dimR = denBK->gCurrentDim( index + 2, NR, TwoSR, IR );
               if ( dimR > 0 ){
                  DimRtotal[ iCenter ] += dimR;
               }
            }
         }
      }
      CenterDims[ iCenter ] = min( DimLtotal[ iCenter ], DimRtotal[ iCenter ] ); // CenterDims contains the min. amount

      //Allocate memory to copy the different parts of the S-object. Use prefactor sqrt((2jR+1)/(2jM+1) * (2jM+1) * (2j+1)) W6J (-1)^(jL+jR+s1+s2) and sum over j.
      if ( CenterDims[ iCenter ] > 0 ){

         // Only if CenterDims[ iCenter ] exists should you allocate the following three arrays
         Lambdas[ iCenter ] = new double[ CenterDims[ iCenter ] ];
              Us[ iCenter ] = new double[ CenterDims[ iCenter ] * DimLtotal[ iCenter ] ];
             VTs[ iCenter ] = new double[ CenterDims[ iCenter ] * DimRtotal[ iCenter ] ];

         const int memsize = DimLtotal[ iCenter ] * DimRtotal[ iCenter ];
         double * mem = new double[ memsize ];
         for ( int cnt = 0; cnt < memsize; cnt++ ){ mem[ cnt ] = 0.0; }

         int dimLtotal2 = 0;
         for ( int NL = SplitSectNM[ iCenter ] - 2; NL <= SplitSectNM[ iCenter ]; NL++ ){
            const int TwoS1 = (( NL + 1 == SplitSectNM[ iCenter ] ) ? 1 : 0 );
            for ( int TwoSL = SplitSectTwoJM[ iCenter ] - TwoS1; TwoSL <= SplitSectTwoJM[ iCenter ] + TwoS1; TwoSL += 2 ){
               if ( TwoSL >= 0 ){
                  const int IL = (( TwoS1 == 1 ) ? Irreps::directProd( Ilocal1, SplitSectIM[ iCenter ] ) : SplitSectIM[ iCenter ] );
                  const int dimL = denBK->gCurrentDim( index, NL, TwoSL, IL );
                  if ( dimL > 0 ){
                     int dimRtotal2 = 0;
                     for ( int NR = SplitSectNM[ iCenter ]; NR <= SplitSectNM[ iCenter ] + 2; NR++ ){
                        const int TwoS2 = (( NR == SplitSectNM[ iCenter ] + 1 ) ? 1 : 0 );
                        for ( int TwoSR = SplitSectTwoJM[ iCenter ] - TwoS2; TwoSR <= SplitSectTwoJM[ iCenter ] + TwoS2; TwoSR += 2 ){
                           if ( TwoSR >= 0 ){
                              const int IR = (( TwoS2 == 1 ) ? Irreps::directProd( Ilocal2, SplitSectIM[ iCenter ] ) : SplitSectIM[ iCenter ] );
                              const int dimR = denBK->gCurrentDim( index + 2, NR, TwoSR, IR );
                              if ( dimR > 0 ){
                                 // Loop over contributing TwoJ's
                                 const int fase = Special::phase( TwoSL + TwoSR + TwoS1 + TwoS2 );
                                 const int TwoJmin = max( abs( TwoSR - TwoSL ), abs( TwoS2 - TwoS1 ) );
                                 const int TwoJmax = min( TwoS1 + TwoS2, TwoSL + TwoSR );
                                 for ( int TwoJ = TwoJmin; TwoJ <= TwoJmax; TwoJ += 2 ){
                                    // Calc prefactor
                                    const double prefactor = fase
                                                           * sqrt( 1.0 * ( TwoJ + 1 ) * ( TwoSR + 1 ) )
                                                           * Wigner::wigner6j( TwoSL, TwoSR, TwoJ, TwoS2, TwoS1, SplitSectTwoJM[ iCenter ] );

                                    // Add them to mem --> += because several TwoJ
                                    double * Block = gStorage( NL, TwoSL, IL, SplitSectNM[ iCenter ] - NL, NR - SplitSectNM[ iCenter ], TwoJ, NR, TwoSR, IR );
                                    for ( int l = 0; l < dimL; l++ ){
                                       for ( int r = 0; r < dimR; r++ ){
                                          mem[ dimLtotal2 + l + DimLtotal[ iCenter ] * ( dimRtotal2 + r ) ] += prefactor * Block[ l + dimL * r ];
                                       }
                                    }
                                 }
                                 dimRtotal2 += dimR;
                              }
                           }
                        }
                     }
                     dimLtotal2 += dimL;
                  }
               }
            }
         }

         // Now mem contains sqrt((2jR+1)/(2jM+1)) * (TT)^{jM nM IM) --> SVD per central symmetry
         char jobz = 'S'; // M x min(M,N) in U and min(M,N) x N in VT
         int lwork = 3 * CenterDims[ iCenter ] + max( max( DimLtotal[ iCenter ], DimRtotal[ iCenter ] ), 4 * CenterDims[ iCenter ] * ( CenterDims[ iCenter ] + 1 ) );
         double * work = new double[ lwork ];
         int * iwork = new int[ 8 * CenterDims[ iCenter ] ];
         int info;

         // dgesdd is not thread-safe in every implementation ( intel MKL is safe, Atlas is not safe )
         #ifndef CHEMPS2_MKL
         #pragma omp critical
         #endif
         dgesdd_( &jobz, DimLtotal + iCenter, DimRtotal + iCenter, mem, DimLtotal + iCenter,
                  Lambdas[ iCenter ], Us[ iCenter ], DimLtotal + iCenter, VTs[ iCenter ], CenterDims + iCenter, work, &lwork, iwork, &info );

         delete [] work;
         delete [] iwork;
         delete [] mem;
      }
   }

   #ifdef CHEMPS2_MPI_COMPILATION
   }
   #endif

   double discardedWeight = 0.0; // Only if change==true; will the discardedWeight be meaningful and different from zero.
   int updateSectors = 0;
   int * NewDims = NULL;

   // If change: determine new virtual dimensions.
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( change ) && ( am_i_master )){
   #else
   if ( change ){
   #endif

      NewDims = new int[ nCenterSectors ];
      // First determine the total number of singular values
      int totalDimSVD = 0;
      for ( int iCenter = 0; iCenter < nCenterSectors; iCenter++ ){
         NewDims[ iCenter ] = CenterDims[ iCenter ];
         totalDimSVD += NewDims[ iCenter ];
      }

      // If larger then the required virtualdimensionD, new virtual dimensions will be set in NewDims.
      if ( totalDimSVD > virtualdimensionD ){
         // Copy them all in 1 array
         double * values = new double[ totalDimSVD ];
         totalDimSVD = 0;
         int inc = 1;
         for ( int iCenter = 0; iCenter < nCenterSectors; iCenter++ ){
            if ( NewDims[ iCenter ] > 0 ){
               dcopy_( NewDims + iCenter, Lambdas[ iCenter ], &inc, values + totalDimSVD, &inc );
               totalDimSVD += NewDims[ iCenter ];
            }
         }

         // Sort them in decreasing order
         char ID = 'D';
         int info;
         dlasrt_( &ID, &totalDimSVD, values, &info ); // Quicksort

         // The D+1'th value becomes the lower bound Schmidt value. Every value smaller than or equal to the D+1'th value is thrown out (hence Dactual <= Ddesired).
         const double lowerBound = values[ virtualdimensionD ];
         for ( int iCenter = 0; iCenter < nCenterSectors; iCenter++ ){
            for ( int cnt = 0; cnt < NewDims[ iCenter ]; cnt++ ){
               if ( Lambdas[ iCenter ][ cnt ] <= lowerBound ){ NewDims[ iCenter ] = cnt; }
            }
         }

         // Discarded weight
         double totalSum = 0.0;
         double discardedSum = 0.0;
         for ( int iCenter = 0; iCenter < nCenterSectors; iCenter++ ){
            for ( int iLocal = 0; iLocal < CenterDims[ iCenter ]; iLocal++ ){
               double temp = ( SplitSectTwoJM[ iCenter ] + 1 ) * Lambdas[ iCenter ][ iLocal ] * Lambdas[ iCenter ][ iLocal ];
               totalSum += temp;
               if ( Lambdas[ iCenter ][ iLocal ] <= lowerBound ){ discardedSum += temp; }
            }
         }
         discardedWeight = discardedSum / totalSum;

         // Clean-up
         delete [] values;
      }

      // Check if there is a sector which differs
      updateSectors = 0;
      for ( int iCenter = 0; iCenter < nCenterSectors; iCenter++ ){
         const int MPSdim = denBK->gCurrentDim( index + 1, SplitSectNM[ iCenter ], SplitSectTwoJM[ iCenter ], SplitSectIM[ iCenter ] );
         if ( NewDims[ iCenter ] != MPSdim ){ updateSectors = 1; }
      }

   }

   #ifdef CHEMPS2_MPI_COMPILATION
   MPIchemps2::broadcast_array_int(    &updateSectors,   1, MPI_CHEMPS2_MASTER );
   MPIchemps2::broadcast_array_double( &discardedWeight, 1, MPI_CHEMPS2_MASTER );
   #endif

   if ( updateSectors == 1 ){

      #ifdef CHEMPS2_MPI_COMPILATION
      if ( NewDims == NULL ){ NewDims = new int[ nCenterSectors ]; }
      MPIchemps2::broadcast_array_int( NewDims, nCenterSectors, MPI_CHEMPS2_MASTER );
      #endif

      for ( int iCenter = 0; iCenter < nCenterSectors; iCenter++ ){
         denBK->SetDim( index + 1, SplitSectNM[ iCenter ], SplitSectTwoJM[ iCenter ], SplitSectIM[ iCenter ], NewDims[ iCenter ] );
      }
      Tleft ->Reset();
      Tright->Reset();

   }

   if ( NewDims != NULL ){ delete [] NewDims; }

   #ifdef CHEMPS2_MPI_COMPILATION
   if ( am_i_master ){
   #endif

   // Copy first dimM per central symmetry sector to the relevant parts
   #pragma omp parallel for schedule(dynamic)
   for ( int iCenter = 0; iCenter < nCenterSectors; iCenter++ ){
      const int dimM = denBK->gCurrentDim( index + 1, SplitSectNM[ iCenter ], SplitSectTwoJM[ iCenter ], SplitSectIM[ iCenter ] );
      if ( dimM > 0 ){
         // U-part: copy
         int dimLtotal2 = 0;
         for ( int NL = SplitSectNM[ iCenter ] - 2; NL <= SplitSectNM[ iCenter ]; NL++ ){
            const int TwoS1 = (( NL + 1 == SplitSectNM[ iCenter ] ) ? 1 : 0 );
            for ( int TwoSL = SplitSectTwoJM[ iCenter ] - TwoS1; TwoSL <= SplitSectTwoJM[ iCenter ] + TwoS1; TwoSL += 2 ){
               if ( TwoSL >= 0 ){
                  const int IL = (( TwoS1 == 1 ) ? Irreps::directProd( Ilocal1, SplitSectIM[ iCenter ] ) : SplitSectIM[ iCenter ] );
                  const int dimL = denBK->gCurrentDim( index, NL, TwoSL, IL );
                  if ( dimL > 0 ){
                     double * TleftBlock = Tleft->gStorage( NL, TwoSL, IL, SplitSectNM[ iCenter ], SplitSectTwoJM[ iCenter ], SplitSectIM[ iCenter ] );
                     const int dimension_limit_right = min( dimM, CenterDims[ iCenter ] );
                     for ( int r = 0; r < dimension_limit_right; r++ ){
                        const double factor = (( movingright ) ? 1.0 : Lambdas[ iCenter ][ r ] );
                        for ( int l = 0; l < dimL; l++ ){
                           TleftBlock[ l + dimL * r ] = factor * Us[ iCenter ][ dimLtotal2 + l + DimLtotal[ iCenter ] * r ];
                        }
                     }
                     for ( int r = dimension_limit_right; r < dimM; r++ ){
                        for ( int l = 0; l < dimL; l++ ){
                           TleftBlock[ l + dimL * r ] = 0.0;
                        }
                     }
                     dimLtotal2 += dimL;
                  }
               }
            }
         }

         // VT-part: copy
         int dimRtotal2 = 0;
         for ( int NR = SplitSectNM[ iCenter ]; NR <= SplitSectNM[ iCenter ] + 2; NR++ ){
            const int TwoS2 = (( NR == SplitSectNM[ iCenter ] + 1 ) ? 1 : 0 );
            for ( int TwoSR = SplitSectTwoJM[ iCenter ] - TwoS2; TwoSR <= SplitSectTwoJM[ iCenter ] + TwoS2; TwoSR += 2 ){
               if ( TwoSR >= 0 ){
                  const int IR = (( TwoS2 == 1 ) ? Irreps::directProd( Ilocal2, SplitSectIM[ iCenter ] ) : SplitSectIM[ iCenter ] );
                  const int dimR = denBK->gCurrentDim( index + 2, NR, TwoSR, IR );
                  if ( dimR > 0 ){
                     double * TrightBlock = Tright->gStorage( SplitSectNM[ iCenter ], SplitSectTwoJM[ iCenter ], SplitSectIM[ iCenter ], NR, TwoSR, IR );
                     const int dimension_limit_left = min( dimM, CenterDims[ iCenter ] );
                     const double factor_base = sqrt( ( SplitSectTwoJM[ iCenter ] + 1.0 ) / ( TwoSR + 1 ) );
                     for ( int l = 0; l < dimension_limit_left; l++ ){
                        const double factor = factor_base * (( movingright ) ? Lambdas[ iCenter ][ l ] : 1.0 );
                        for ( int r = 0; r < dimR; r++ ){
                           TrightBlock[ l + dimM * r ] = factor * VTs[ iCenter ][ l + CenterDims[ iCenter ] * ( dimRtotal2 + r ) ];
                        }
                     }
                     for ( int r = 0; r < dimR; r++ ){
                        for ( int l = dimension_limit_left; l < dimM; l++ ){
                           TrightBlock[ l + dimM * r ] = 0.0;
                        }
                     }
                     dimRtotal2 += dimR;
                  }
               }
            }
         }
      }
   }

   #ifdef CHEMPS2_MPI_COMPILATION
   }
   MPIchemps2::broadcast_tensor( Tleft,  MPI_CHEMPS2_MASTER );
   MPIchemps2::broadcast_tensor( Tright, MPI_CHEMPS2_MASTER );
   #endif

   // Clean up
   delete [] SplitSectNM;
   delete [] SplitSectTwoJM;
   delete [] SplitSectIM;
   #ifdef CHEMPS2_MPI_COMPILATION
   if ( am_i_master )
   #endif
   {
      for ( int iCenter = 0; iCenter < nCenterSectors; iCenter++ ){
         if ( CenterDims[ iCenter ] > 0 ){
            delete []      Us[ iCenter ];
            delete [] Lambdas[ iCenter ];
            delete []     VTs[ iCenter ];
         }
      }
      delete [] Us;
      delete [] Lambdas;
      delete [] VTs;
      delete [] CenterDims;
      delete [] DimLtotal;
      delete [] DimRtotal;
   }

   return discardedWeight;

}

void CheMPS2::Sobject::prog2symm(){

   #pragma omp parallel for schedule(dynamic)
   for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){

      int dim = kappa2index[ ikappa + 1 ] - kappa2index[ ikappa ];
      double alpha = sqrt( sectorTwoSR[ ikappa ] + 1.0 );
      int inc = 1;
      dscal_( &dim, &alpha, storage + kappa2index[ ikappa ], &inc );

   }

}

void CheMPS2::Sobject::symm2prog(){

   #pragma omp parallel for schedule(dynamic)
   for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){

      int dim = kappa2index[ ikappa + 1 ] - kappa2index[ ikappa ];
      double alpha = 1.0 / sqrt( sectorTwoSR[ ikappa ] + 1.0 );
      int inc = 1;
      dscal_( &dim, &alpha, storage + kappa2index[ ikappa ], &inc );

   }

}

void CheMPS2::Sobject::addNoise( const double NoiseLevel ){
   
   for ( int cnt = 0; cnt < gKappa2index( gNKappa() ); cnt++ ){
      const double RN = ( ( double ) rand() ) / RAND_MAX - 0.5;
      gStorage()[ cnt ] += RN * NoiseLevel;
   }

}


