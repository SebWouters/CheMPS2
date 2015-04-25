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

#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include "DMRGSCFVmatRotations.h"
#include "Lapack.h"

using std::min;
using std::max;

CheMPS2::DMRGSCFVmatRotations::DMRGSCFVmatRotations(Hamiltonian * HamOrigIn, DMRGSCFindices * iHandlerIn){

   HamOrig = HamOrigIn;
   iHandler = iHandlerIn;
   SymmInfo.setGroup(HamOrig->getNGroup());
   numberOfIrreps = SymmInfo.getNumberOfIrreps();

}

CheMPS2::DMRGSCFVmatRotations::~DMRGSCFVmatRotations(){ }

void CheMPS2::DMRGSCFVmatRotations::fillVmatRotated(FourIndex * VmatRotated, DMRGSCFunitary * unitary, double * temp1, double * temp2) const{
   
   //Two-body terms --> use eightfold permutation symmetry
   for (int irrep1 = 0; irrep1<numberOfIrreps; irrep1++){
      for (int irrep2 = irrep1; irrep2<numberOfIrreps; irrep2++){
         const int productSymm = Irreps::directProd(irrep1,irrep2);
         for (int irrep3 = irrep1; irrep3<numberOfIrreps; irrep3++){
            const int irrep4 = Irreps::directProd(productSymm,irrep3);
            if (irrep4>=irrep2){
            
               int linsize1 = iHandler->getNORB(irrep1);
               int linsize2 = iHandler->getNORB(irrep2);
               int linsize3 = iHandler->getNORB(irrep3);
               int linsize4 = iHandler->getNORB(irrep4);
               
               if ((linsize1>0) && (linsize2>0) && (linsize3>0) && (linsize4>0)){
                  
                  for (int cnt1=0; cnt1<linsize1; cnt1++){
                     for (int cnt2=0; cnt2<linsize2; cnt2++){
                        for (int cnt3=0; cnt3<linsize3; cnt3++){
                           for (int cnt4=0; cnt4<linsize4; cnt4++){
                              temp1[cnt1 + linsize1 * ( cnt2 + linsize2 * (cnt3 + linsize3 * cnt4) ) ]
                                = HamOrig->getVmat( iHandler->getOrigNOCCstart(irrep1) + cnt1, iHandler->getOrigNOCCstart(irrep2) + cnt2,
                                                    iHandler->getOrigNOCCstart(irrep3) + cnt3, iHandler->getOrigNOCCstart(irrep4) + cnt4 );
                           }
                        }
                     }
                  }
                  
                  char trans = 'T';
                  char notra = 'N';
                  double alpha = 1.0;
                  double beta = 0.0; //SET !!!
                  
                  int rightdim = linsize2 * linsize3 * linsize4; //(ijkl) -> (ajkl)
                  dgemm_(&notra,&notra,&linsize1,&rightdim,&linsize1, &alpha, unitary->getBlock(irrep1),&linsize1,temp1,&linsize1, &beta,temp2,&linsize1);
                  
                  int leftdim = linsize1 * linsize2 * linsize3; //(ajkl) -> (ajkd)
                  dgemm_(&notra,&trans,&leftdim,&linsize4,&linsize4, &alpha,temp2,&leftdim, unitary->getBlock(irrep4),&linsize4, &beta,temp1,&leftdim);
                  
                  int jump = leftdim; //(ajkd) -> (ajcd)
                  leftdim = linsize1 * linsize2;
                  for (int bla=0; bla<linsize4; bla++){
                     dgemm_(&notra,&trans,&leftdim,&linsize3,&linsize3, &alpha,temp1+jump*bla,&leftdim, unitary->getBlock(irrep3),&linsize3, &beta,temp2+jump*bla,&leftdim);
                  }
                  
                  jump = leftdim;
                  rightdim = linsize3*linsize4;
                  for (int bla=0; bla<rightdim; bla++){
                     dgemm_(&notra,&trans,&linsize1,&linsize2,&linsize2,&alpha,temp2+jump*bla,&linsize1, unitary->getBlock(irrep2),&linsize2,&beta,temp1+jump*bla,&linsize1);
                  }
                  
                  for (int cnt1=0; cnt1<linsize1; cnt1++){
                     for (int cnt2=0; cnt2<linsize2; cnt2++){
                        for (int cnt3=0; cnt3<linsize3; cnt3++){
                           for (int cnt4=0; cnt4<linsize4; cnt4++){
                              VmatRotated->set(irrep1, irrep2, irrep3, irrep4, cnt1, cnt2, cnt3, cnt4,
                                               temp1[cnt1 + linsize1 * ( cnt2 + linsize2 * (cnt3 + linsize3 * cnt4) ) ] );
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

void CheMPS2::DMRGSCFVmatRotations::fillVmatDMRG(Hamiltonian * HamDMRG, DMRGSCFunitary * unitary, double * temp1, double * temp2) const{
   
   //Two-body terms --> use eightfold permutation symmetry in the irreps :-)
   for (int irrep1 = 0; irrep1<numberOfIrreps; irrep1++){
      for (int irrep2 = irrep1; irrep2<numberOfIrreps; irrep2++){
         const int productSymm = Irreps::directProd(irrep1,irrep2);
         for (int irrep3 = irrep1; irrep3<numberOfIrreps; irrep3++){
            const int irrep4 = Irreps::directProd(productSymm,irrep3);
            if (irrep4>=irrep2){
            
               int linsizeDMRG1 = iHandler->getNDMRG(irrep1);
               int linsizeDMRG2 = iHandler->getNDMRG(irrep2);
               int linsizeDMRG3 = iHandler->getNDMRG(irrep3);
               int linsizeDMRG4 = iHandler->getNDMRG(irrep4);
               
               if ((linsizeDMRG1>0) && (linsizeDMRG2>0) && (linsizeDMRG3>0) && (linsizeDMRG4>0)){
               
                  int linsizeORIG1 = iHandler->getNORB(irrep1);
                  int linsizeORIG2 = iHandler->getNORB(irrep2);
                  int linsizeORIG3 = iHandler->getNORB(irrep3);
                  int linsizeORIG4 = iHandler->getNORB(irrep4);
                  
                  for (int cnt1=0; cnt1<linsizeORIG1; cnt1++){
                     for (int cnt2=0; cnt2<linsizeORIG2; cnt2++){
                        for (int cnt3=0; cnt3<linsizeORIG3; cnt3++){
                           for (int cnt4=0; cnt4<linsizeORIG4; cnt4++){
                              temp1[cnt1 + linsizeORIG1 * ( cnt2 + linsizeORIG2 * (cnt3 + linsizeORIG3 * cnt4) ) ]
                                = HamOrig->getVmat( iHandler->getOrigNOCCstart(irrep1) + cnt1, iHandler->getOrigNOCCstart(irrep2) + cnt2,
                                                    iHandler->getOrigNOCCstart(irrep3) + cnt3, iHandler->getOrigNOCCstart(irrep4) + cnt4 );
                           }
                        }
                     }
                  }
                  
                  char trans = 'T';
                  char notra = 'N';
                  double alpha = 1.0;
                  double beta  = 0.0; //SET !!!
                  
                  int rightdim = linsizeORIG2 * linsizeORIG3 * linsizeORIG4; //(ijkl) -> (ajkl)
                  double * Umx = unitary->getBlock(irrep1) + iHandler->getNOCC(irrep1);
                  dgemm_(&notra, &notra, &linsizeDMRG1, &rightdim, &linsizeORIG1, &alpha, Umx, &linsizeORIG1, temp1, &linsizeORIG1, &beta, temp2, &linsizeDMRG1);
                  
                  int leftdim = linsizeDMRG1 * linsizeORIG2 * linsizeORIG3; //(ajkl) -> (ajkd)
                  Umx = unitary->getBlock(irrep4) + iHandler->getNOCC(irrep4);
                  dgemm_(&notra, &trans, &leftdim, &linsizeDMRG4, &linsizeORIG4, &alpha, temp2, &leftdim, Umx, &linsizeORIG4, &beta, temp1, &leftdim);
                  
                  int jump1 = linsizeDMRG1 * linsizeORIG2 * linsizeORIG3; //(ajkd) -> (ajcd)
                  int jump2 = linsizeDMRG1 * linsizeORIG2 * linsizeDMRG3;
                  leftdim   = linsizeDMRG1 * linsizeORIG2;
                  Umx = unitary->getBlock(irrep3) + iHandler->getNOCC(irrep3);
                  for (int bla=0; bla<linsizeDMRG4; bla++){
                     dgemm_(&notra, &trans, &leftdim, &linsizeDMRG3, &linsizeORIG3, &alpha, temp1+jump1*bla, &leftdim, Umx, &linsizeORIG3, &beta, temp2+jump2*bla, &leftdim);
                  }
                  
                  jump2    = linsizeDMRG1 * linsizeORIG2;
                  jump1    = linsizeDMRG1 * linsizeDMRG2;
                  rightdim = linsizeDMRG3 * linsizeDMRG4;
                  Umx = unitary->getBlock(irrep2) + iHandler->getNOCC(irrep2);
                  for (int bla=0; bla<rightdim; bla++){
                     dgemm_(&notra, &trans, &linsizeDMRG1, &linsizeDMRG2, &linsizeORIG2, &alpha, temp2+jump2*bla, &linsizeDMRG1, Umx, &linsizeORIG2, &beta, temp1+jump1*bla, &linsizeDMRG1);
                  }
                  
                  for (int cnt1=0; cnt1<linsizeDMRG1; cnt1++){
                     for (int cnt2=0; cnt2<linsizeDMRG2; cnt2++){
                        for (int cnt3=0; cnt3<linsizeDMRG3; cnt3++){
                           for (int cnt4=0; cnt4<linsizeDMRG4; cnt4++){
                              HamDMRG->setVmat( iHandler->getDMRGcumulative(irrep1) + cnt1, iHandler->getDMRGcumulative(irrep2) + cnt2,
                                                iHandler->getDMRGcumulative(irrep3) + cnt3, iHandler->getDMRGcumulative(irrep4) + cnt4,
                                                temp1[cnt1 + linsizeDMRG1 * ( cnt2 + linsizeDMRG2 * (cnt3 + linsizeDMRG3 * cnt4) ) ] );
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

void CheMPS2::DMRGSCFVmatRotations::fillRotatedTEI(DMRGSCFintegrals * theRotatedTEI, DMRGSCFunitary * unitary, double * temp1, double * temp2) const{

   // First do Coulomb object : ( c1 <= c2 | a1 <= a2 )
   for (int Ic1 = 0; Ic1 < numberOfIrreps; Ic1++){
      for (int Ic2 = Ic1; Ic2 < numberOfIrreps; Ic2++){
         const int Icc = Irreps::directProd( Ic1, Ic2 );
         for (int Ia1 = 0; Ia1 < numberOfIrreps; Ia1++){
            const int Ia2 = Irreps::directProd( Ia1, Icc );
            if ( Ia1 <= Ia2 ){
               
               int linsize_orig_c1 = iHandler->getNORB( Ic1 );
               int linsize_orig_c2 = iHandler->getNORB( Ic2 );
               int linsize_small_c1 = iHandler->getNOCC( Ic1 ) + iHandler->getNDMRG( Ic1 );
               int linsize_small_c2 = iHandler->getNOCC( Ic2 ) + iHandler->getNDMRG( Ic2 );
               int linsize_a1 = iHandler->getNORB( Ia1 );
               int linsize_a2 = iHandler->getNORB( Ia2 );
               
               if (( linsize_small_c1 > 0 ) && ( linsize_small_c2 > 0 ) && ( linsize_a1 > 0 ) && ( linsize_a2 > 0 )){
               
                  for (int c1 = 0; c1 < linsize_orig_c1; c1++){
                     for (int c2 = 0; c2 < linsize_orig_c2; c2++){
                        for (int a1 = 0; a1 < linsize_a1; a1++){
                           for (int a2 = 0; a2 < linsize_a2; a2++){
                              // We try to make the Coulomb elements !
                              temp1[ c1 + linsize_orig_c1 * ( c2 + linsize_orig_c2 * ( a1 + linsize_a1 * a2 ) ) ]
                                = HamOrig->getVmat( iHandler->getOrigNOCCstart( Ic1 ) + c1, iHandler->getOrigNOCCstart( Ia1 ) + a1,
                                                    iHandler->getOrigNOCCstart( Ic2 ) + c2, iHandler->getOrigNOCCstart( Ia2 ) + a2 );
                           }
                        }
                     }
                  }
                  
                  char trans = 'T';
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0; //SET !!!
                  
                  int rightdim = linsize_orig_c2 * linsize_a1 * linsize_a2; //( ij | kl ) -> ( c1 j | kl )
                  double * Umx = unitary->getBlock( Ic1 );
                  dgemm_(&notrans, &notrans, &linsize_small_c1, &rightdim, &linsize_orig_c1, &alpha,
                         Umx, &linsize_orig_c1, temp1, &linsize_orig_c1, &beta, temp2, &linsize_small_c1);
                  
                  rightdim = linsize_a1 * linsize_a2; //( c1 j | kl ) -> ( c1 c2 | kl )
                  int jump1 = linsize_small_c1 * linsize_small_c2;
                  int jump2 = linsize_small_c1 * linsize_orig_c2;
                  Umx = unitary->getBlock( Ic2 );
                  for (int loop = 0; loop < rightdim; loop++){
                     dgemm_(&notrans, &trans, &linsize_small_c1, &linsize_small_c2, &linsize_orig_c2, &alpha,
                            temp2+loop*jump2, &linsize_small_c1, Umx, &linsize_orig_c2, &beta, temp1+loop*jump1, &linsize_small_c1);
                  }
                  
                  int leftdim = linsize_small_c1 * linsize_small_c2 * linsize_a1; //( c1 c2 | kl ) -> ( c1 c2 | k a2 )
                  Umx = unitary->getBlock( Ia2 );
                  dgemm_(&notrans, &trans, &leftdim, &linsize_a2, &linsize_a2, &alpha,
                         temp1, &leftdim, Umx, &linsize_a2, &beta, temp2, &leftdim);
                  
                  jump1 = linsize_small_c1 * linsize_small_c2 * linsize_a1; //( c1 c2 | k a2 ) -> ( c1 c2 | a1 a2 )
                  leftdim = linsize_small_c1 * linsize_small_c2;
                  Umx = unitary->getBlock( Ia1 );
                  for (int loop = 0; loop < linsize_a2; loop++){
                     dgemm_(&notrans, &trans, &leftdim, &linsize_a1, &linsize_a1, &alpha,
                            temp2+jump1*loop, &leftdim, Umx, &linsize_a1, &beta, temp1+jump1*loop, &leftdim);
                  }
                  
                  for (int c1 = 0; c1 < linsize_small_c1; c1++){
                     for (int c2 = 0; c2 < linsize_small_c2; c2++){
                        for (int a1 = 0; a1 < linsize_a1; a1++){
                           for (int a2 = 0; a2 < linsize_a2; a2++){
                              theRotatedTEI->set_coulomb( Ic1, Ic2, Ia1, Ia2, c1, c2, a1, a2,
                                  temp1[ c1 + linsize_small_c1 * ( c2 + linsize_small_c2 * ( a1 + linsize_a1 * a2 ) ) ] );
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   // Now do Exchange object ( c1 v1 | c2 v2 ) with c1 <= c2
   for (int Ic1 = 0; Ic1 < numberOfIrreps; Ic1++){
      for (int Ic2 = Ic1; Ic2 < numberOfIrreps; Ic2++){
         const int Icc = Irreps::directProd( Ic1, Ic2 );
         for (int Iv1 = 0; Iv1 < numberOfIrreps; Iv1++){
            const int Iv2 = Irreps::directProd( Iv1, Icc );
               
            int linsize_orig_c1 = iHandler->getNORB( Ic1 );
            int linsize_orig_c2 = iHandler->getNORB( Ic2 );
            int linsize_small_c1 = iHandler->getNOCC( Ic1 ) + iHandler->getNDMRG( Ic1 );
            int linsize_small_c2 = iHandler->getNOCC( Ic2 ) + iHandler->getNDMRG( Ic2 );
            
            int linsize_orig_v1 = iHandler->getNORB( Iv1 );
            int linsize_orig_v2 = iHandler->getNORB( Iv2 );
            int linsize_small_v1 = iHandler->getNVIRT( Iv1 );
            int linsize_small_v2 = iHandler->getNVIRT( Iv2 );
            
            if (( linsize_small_c1 > 0 ) && ( linsize_small_c2 > 0 ) && ( linsize_small_v1 > 0 ) && ( linsize_small_v2 > 0 )){
            
               for (int c1 = 0; c1 < linsize_orig_c1; c1++){
                  for (int c2 = 0; c2 < linsize_orig_c2; c2++){
                     for (int v1 = 0; v1 < linsize_orig_v1; v1++){
                        for (int v2 = 0; v2 < linsize_orig_v2; v2++){
                           // We try to make the Exchange elements !
                           temp1[ c1 + linsize_orig_c1 * ( c2 + linsize_orig_c2 * ( v1 + linsize_orig_v1 * v2 ) ) ]
                             = HamOrig->getVmat( iHandler->getOrigNOCCstart( Ic1 ) + c1, iHandler->getOrigNOCCstart( Ic2 ) + c2,
                                                 iHandler->getOrigNOCCstart( Iv1 ) + v1, iHandler->getOrigNOCCstart( Iv2 ) + v2 );
                        }
                     }
                  }
               }
               
               char trans = 'T';
               char notrans = 'N';
               double alpha = 1.0;
               double beta = 0.0; //SET !!!
               const int shiftv1 = linsize_orig_v1 - linsize_small_v1;
               const int shiftv2 = linsize_orig_v2 - linsize_small_v2;
               
               int rightdim = linsize_orig_c2 * linsize_orig_v1 * linsize_orig_v2; //( ij | kl ) -> ( c1 j | kl )
               double * Umx = unitary->getBlock( Ic1 );
               dgemm_(&notrans, &notrans, &linsize_small_c1, &rightdim, &linsize_orig_c1, &alpha,
                      Umx, &linsize_orig_c1, temp1, &linsize_orig_c1, &beta, temp2, &linsize_small_c1);
               
               rightdim = linsize_orig_v1 * linsize_orig_v2; //( c1 j | kl ) -> ( c1 j | c2 l )
               int jump1 = linsize_small_c1 * linsize_small_c2;
               int jump2 = linsize_small_c1 * linsize_orig_c2;
               Umx = unitary->getBlock( Ic2 );
               for (int loop = 0; loop < rightdim; loop++){
                  dgemm_(&notrans, &trans, &linsize_small_c1, &linsize_small_c2, &linsize_orig_c2, &alpha,
                         temp2+loop*jump2, &linsize_small_c1, Umx, &linsize_orig_c2, &beta, temp1+loop*jump1, &linsize_small_c1);
               }
               
               int leftdim = linsize_small_c1 * linsize_small_c2 * linsize_orig_v1; //( c1 j | c2 l ) -> ( c1 j | c2 v2 )
               Umx = unitary->getBlock( Iv2 ) + shiftv2;
               dgemm_(&notrans, &trans, &leftdim, &linsize_small_v2, &linsize_orig_v2, &alpha,
                      temp1, &leftdim, Umx, &linsize_orig_v2, &beta, temp2, &leftdim);
               
               jump1 = linsize_small_c1 * linsize_small_c2 * linsize_orig_v1; //( c1 j | c2 v2 ) ->  ( c1 v1 | c2 v2 )
               jump2 = linsize_small_c1 * linsize_small_c2 * linsize_small_v1;
               leftdim = linsize_small_c1 * linsize_small_c2;
               Umx = unitary->getBlock( Iv1 ) + shiftv1;
               for (int loop = 0; loop < linsize_small_v2; loop++){
                  dgemm_(&notrans, &trans, &leftdim, &linsize_small_v1, &linsize_orig_v1, &alpha,
                         temp2+jump1*loop, &leftdim, Umx, &linsize_orig_v1, &beta, temp1+jump2*loop, &leftdim);
               }
               
               for (int c1 = 0; c1 < linsize_small_c1; c1++){
                  for (int c2 = 0; c2 < linsize_small_c2; c2++){
                     for (int v1 = 0; v1 < linsize_small_v1; v1++){
                        for (int v2 = 0; v2 < linsize_small_v2; v2++){
                           theRotatedTEI->set_exchange( Ic1, Ic2, Iv1, Iv2, c1, c2, shiftv1 + v1, shiftv2 + v2,
                               temp1[ c1 + linsize_small_c1 * ( c2 + linsize_small_c2 * ( v1 + linsize_small_v1 * v2 ) ) ] );
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::DMRGSCFVmatRotations::fillVmatRotatedBlockWise(FourIndex * VmatRotated, DMRGSCFunitary * unitary, double * mem1, double * mem2, double * mem3, const int maxBlockSize, const bool cutCorners) const{
   
   /***********************************************************
   **   Two-body terms; use eightfold permutation symmetry   **
   **   Requires 3 buffers of size maxBlockSize^4            **
   **   maxBlockSize = ceil(maxlinsize/factor) at alloc      **
   ***********************************************************/
   for (int irrep1 = 0; irrep1<numberOfIrreps; irrep1++){
      for (int irrep2 = irrep1; irrep2<numberOfIrreps; irrep2++){
         const int productSymm = Irreps::directProd(irrep1,irrep2);
         for (int irrep3 = irrep1; irrep3<numberOfIrreps; irrep3++){
            const int irrep4 = Irreps::directProd(productSymm,irrep3);
            if (irrep4>=irrep2){
            
               const int linsize1 = iHandler->getNORB(irrep1);
               const int linsize2 = iHandler->getNORB(irrep2);
               const int linsize3 = iHandler->getNORB(irrep3);
               const int linsize4 = iHandler->getNORB(irrep4);
               
               const bool caseA = ((irrep1 == irrep2) && (irrep1 == irrep3) && (irrep1 == irrep4)); //I1 = I2 = I3 = I4
               const bool caseB = ((!caseA) && (irrep1 == irrep2) && (irrep3 == irrep4));           //I1 = I2 < I3 = I4
               const bool caseC = ((!caseA) && (irrep1 == irrep3) && (irrep2 == irrep4));           //I1 = I3 < I2 = I4
               //caseD : All irreps different
               /* if (cutCorners):
                     For the Gradient and Hessian, not all matrix elements are required: Maximum two indices are virtual.
                     Thereto split up all target indices in an OA (occupied + active) part and a V (virtual) part and see that max. two indices are virtual.

                     CASE A: if I1 = I2 = I3 = I4 the following combinations are allowed
                        0: OA  OA  OA  OA
                        1: OA  OA  OA  V
                        1: OA  OA  V   OA
                        2: OA  OA  V   V
                        2: OA  V   OA  V
                     due to either ( i1 < i2 <= i4 and i1 <= i3 ) or ( i1 = i2 <= i3 <= i4 ) .

                     CASE B: if I1 = I2 < I3 = I4 the following combinations are allowed
                        0: OA  OA  OA  OA
                        1: OA  OA  OA  V
                        1: OA  OA  V   OA
                        1: OA  V   OA  OA
                        2: OA  OA  V   V
                        2: OA  V   V   OA
                        2: OA  V   OA  V
                        2: V   V   OA  OA
                     due to either ( i1 < i2 and i3,i4 all ) or ( i1 = i2 <= i3 <= i4 ) .

                     CASE C: if I1 = I3 < I2 = I4 the following combinations are allowed
                        0: OA  OA  OA  OA
                        1: OA  OA  OA  V
                        1: OA  OA  V   OA
                        2: OA  OA  V   V
                        2: V   OA  V   OA
                        2: OA  V   OA  V
                     due to either ( i1 <= i3 and i2 <= i4 ) .

                     CASE D: all irreps different: max two blocks virtual .
               */
               
               if ((linsize1>0) && (linsize2>0) && (linsize3>0) && (linsize4>0)){
               
                  /** The scheme: 
                      Split up each of the linsizes into factor parts of size (linsize/factor). Make sure the boundaries are clear.
                      
                      Algorithm:
                         Reset the FourIndex object for the particular symmetry case under consideration --> Vmat(i,j,k,l) = 0.0
                         Loop block1, block2, block3, block4 --> factor^4:
                            Loop block1 --> factor:
                               Copy HamOrig->getVmat(i,j,k,l) for the particular block into mem2_ijkl --> cost (linsize/factor)^4
                               Rotate mem1_ajkl = U_ai mem2_ijkl (but only block1 indices) --> cost(linsize/factor)^5
                               Loop block2 --> factor:
                                  Rotate mem2_abkl = U_bj mem1_ajkl (only block2 indices) --> cost (linsize/factor)^5
                                  Loop block3 --> factor:
                                     Rotate mem3_abcl = U_ck mem2_abkl (only block3 indices) --> cost (linsize/factor)^5
                                     Loop block4 --> factor:
                                        Calculate V_partial_abcd = U_dl mem3_abcl (only block4 indices) --> cost (linsize/factor)^5
                                        Add (not set!!!) V_partial_abcd to FourIndex (for the different ijkl-blocks) --> cost (linsize/factor)^4
                      Total cost: f^4 * ( (l/f)^4 + f * ( (l/f)^5 + f * ( (l/f)^5 + f * ( (l/f)^5 + f * ( (l/f)^5 + (l/f)^4 ) ) ) ) ) = f^3 * l^5
                  */
               
                  //Split up the linsizes for the original Hamiltonian
                  int factor1 = max( (int) ( ceil((1.0 * linsize1) / maxBlockSize) + 0.01 ) , 1 ); //factor >= linsize/maxBlockSize
                  int factor2 = max( (int) ( ceil((1.0 * linsize2) / maxBlockSize) + 0.01 ) , 1 );
                  int factor3 = max( (int) ( ceil((1.0 * linsize3) / maxBlockSize) + 0.01 ) , 1 );
                  int factor4 = max( (int) ( ceil((1.0 * linsize4) / maxBlockSize) + 0.01 ) , 1 );
                  
                  const int blocksize1 = min( (int) ( ceil( (1.0 * linsize1) / factor1 ) + 0.01 ) , maxBlockSize ); //Hence at most maxBlockSize
                  const int blocksize2 = min( (int) ( ceil( (1.0 * linsize2) / factor2 ) + 0.01 ) , maxBlockSize );
                  const int blocksize3 = min( (int) ( ceil( (1.0 * linsize3) / factor3 ) + 0.01 ) , maxBlockSize );
                  const int blocksize4 = min( (int) ( ceil( (1.0 * linsize4) / factor4 ) + 0.01 ) , maxBlockSize );
                  
                  while (factor1 * blocksize1 < linsize1){ factor1++; }
                  while (factor2 * blocksize2 < linsize2){ factor2++; }
                  while (factor3 * blocksize3 < linsize3){ factor3++; }
                  while (factor4 * blocksize4 < linsize4){ factor4++; }
               
                  //Split up the linsizes for the rotated Hamiltonian
                  const int linsizeV1  = iHandler->getNVIRT(irrep1);
                  const int linsizeOA1 = linsize1 - linsizeV1;
                  const int linsizeV2  = iHandler->getNVIRT(irrep2);
                  const int linsizeOA2 = linsize2 - linsizeV2;
                  const int linsizeV3  = iHandler->getNVIRT(irrep3);
                  const int linsizeOA3 = linsize3 - linsizeV3;
                  const int linsizeV4  = iHandler->getNVIRT(irrep4);
                  const int linsizeOA4 = linsize4 - linsizeV4;
                  
                  int factorV1 = (linsizeV1 == 0) ? 0 : max( (int) ( ceil((1.0 * linsizeV1) / maxBlockSize) + 0.01 ) , 1 );
                  int factorV2 = (linsizeV2 == 0) ? 0 : max( (int) ( ceil((1.0 * linsizeV2) / maxBlockSize) + 0.01 ) , 1 );
                  int factorV3 = (linsizeV3 == 0) ? 0 : max( (int) ( ceil((1.0 * linsizeV3) / maxBlockSize) + 0.01 ) , 1 );
                  int factorV4 = (linsizeV4 == 0) ? 0 : max( (int) ( ceil((1.0 * linsizeV4) / maxBlockSize) + 0.01 ) , 1 );
                  
                  const int blocksizeV1 = (linsizeV1 == 0) ? 1 : min( (int) ( ceil( (1.0 * linsizeV1) / factorV1 ) + 0.01 ) , maxBlockSize );
                  const int blocksizeV2 = (linsizeV2 == 0) ? 1 : min( (int) ( ceil( (1.0 * linsizeV2) / factorV2 ) + 0.01 ) , maxBlockSize );
                  const int blocksizeV3 = (linsizeV3 == 0) ? 1 : min( (int) ( ceil( (1.0 * linsizeV3) / factorV3 ) + 0.01 ) , maxBlockSize );
                  const int blocksizeV4 = (linsizeV4 == 0) ? 1 : min( (int) ( ceil( (1.0 * linsizeV4) / factorV4 ) + 0.01 ) , maxBlockSize );
                  
                  if (linsizeV1 > 0){ while (factorV1 * blocksizeV1 < linsizeV1){ factorV1++; } }
                  if (linsizeV2 > 0){ while (factorV2 * blocksizeV2 < linsizeV2){ factorV2++; } }
                  if (linsizeV3 > 0){ while (factorV3 * blocksizeV3 < linsizeV3){ factorV3++; } }
                  if (linsizeV4 > 0){ while (factorV4 * blocksizeV4 < linsizeV4){ factorV4++; } }
                  
                  int factorOA1 = (linsizeOA1 == 0) ? 0 : max( (int) ( ceil((1.0 * linsizeOA1) / maxBlockSize) + 0.01 ) , 1 );
                  int factorOA2 = (linsizeOA2 == 0) ? 0 : max( (int) ( ceil((1.0 * linsizeOA2) / maxBlockSize) + 0.01 ) , 1 );
                  int factorOA3 = (linsizeOA3 == 0) ? 0 : max( (int) ( ceil((1.0 * linsizeOA3) / maxBlockSize) + 0.01 ) , 1 );
                  int factorOA4 = (linsizeOA4 == 0) ? 0 : max( (int) ( ceil((1.0 * linsizeOA4) / maxBlockSize) + 0.01 ) , 1 );
                  
                  const int blocksizeOA1 = (linsizeOA1 == 0) ? 1 : min( (int) ( ceil( (1.0 * linsizeOA1) / factorOA1 ) + 0.01 ) , maxBlockSize );
                  const int blocksizeOA2 = (linsizeOA2 == 0) ? 1 : min( (int) ( ceil( (1.0 * linsizeOA2) / factorOA2 ) + 0.01 ) , maxBlockSize );
                  const int blocksizeOA3 = (linsizeOA3 == 0) ? 1 : min( (int) ( ceil( (1.0 * linsizeOA3) / factorOA3 ) + 0.01 ) , maxBlockSize );
                  const int blocksizeOA4 = (linsizeOA4 == 0) ? 1 : min( (int) ( ceil( (1.0 * linsizeOA4) / factorOA4 ) + 0.01 ) , maxBlockSize );
                  
                  if (linsizeOA1 > 0){ while (factorOA1 * blocksizeOA1 < linsizeOA1){ factorOA1++; } }
                  if (linsizeOA2 > 0){ while (factorOA2 * blocksizeOA2 < linsizeOA2){ factorOA2++; } }
                  if (linsizeOA3 > 0){ while (factorOA3 * blocksizeOA3 < linsizeOA3){ factorOA3++; } }
                  if (linsizeOA4 > 0){ while (factorOA4 * blocksizeOA4 < linsizeOA4){ factorOA4++; } }
                  
                  //Reset the FourIndex object
                  for (int cnt1=0; cnt1<linsize1; cnt1++){
                     for (int cnt2=0; cnt2<linsize2; cnt2++){
                        for (int cnt3=0; cnt3<linsize3; cnt3++){
                           for (int cnt4=0; cnt4<linsize4; cnt4++){
                              VmatRotated->set(irrep1, irrep2, irrep3, irrep4, cnt1, cnt2, cnt3, cnt4, 0.0);
                           }
                        }
                     }
                  }
                  
                  //Loop original blocks
                  for (int origBlock1=0; origBlock1<factor1; origBlock1++){
                  
                     const int origStart1 = origBlock1 * blocksize1;
                     const int origStop1  = min( (origBlock1 + 1) * blocksize1 , linsize1 );
                     const int origSize1  = max( origStop1 - origStart1 , 0 );
                     
                     for (int origBlock2=0; origBlock2<factor2; origBlock2++){
                     
                        const int origStart2 = origBlock2 * blocksize2;
                        const int origStop2  = min( (origBlock2 + 1) * blocksize2 , linsize2 );
                        const int origSize2  = max( origStop2 - origStart2 , 0 );
                     
                        for (int origBlock3=0; origBlock3<factor3; origBlock3++){
                        
                           const int origStart3 = origBlock3 * blocksize3;
                           const int origStop3  = min( (origBlock3 + 1) * blocksize3 , linsize3 );
                           const int origSize3  = max( origStop3 - origStart3 , 0 );
                        
                           for (int origBlock4=0; origBlock4<factor4; origBlock4++){
                           
                              const int origStart4 = origBlock4 * blocksize4;
                              const int origStop4  = min( (origBlock4 + 1) * blocksize4 , linsize4 );
                              const int origSize4  = max( origStop4 - origStart4 , 0 );
                              
                              //Loop the target block 1
                              int startTB1 = 0;
                              int stopTB1  = factorOA1 + factorV1;
                              if ((cutCorners) && (caseA)){ stopTB1 = factorOA1; } //Only for case A the first index CANNOT be V.
                              for (int targetBlock1 = startTB1; targetBlock1 < stopTB1; targetBlock1++){
                  
                                 const bool block1_OA   = (targetBlock1 < factorOA1) ? true : false;
                                 const int targetStart1 = (block1_OA) ? (targetBlock1 * blocksizeOA1)
                                                                      : (linsizeOA1 + (targetBlock1 - factorOA1) * blocksizeV1);
                                 const int targetStop1  = (block1_OA) ? min( targetStart1 + blocksizeOA1 , linsizeOA1 )
                                                                      : min( targetStart1 + blocksizeV1  , linsize1   );
                                 const int targetSize1  = max( targetStop1 - targetStart1 , 0 );
                                 
                                 //Copy HamOrig->getVmat(i,j,k,l) for the particular ORIGINAL block into mem2_ijkl
                                 for (int origIndex1=0; origIndex1<origSize1; origIndex1++){
                                    for (int origIndex2=0; origIndex2<origSize2; origIndex2++){
                                       for (int origIndex3=0; origIndex3<origSize3; origIndex3++){
                                          for (int origIndex4=0; origIndex4<origSize4; origIndex4++){
                                             mem2[ origIndex1 + origSize1 * (origIndex2 + origSize2 * (origIndex3 + origSize3 * origIndex4) ) ]
                                                = HamOrig->getVmat( iHandler->getOrigNOCCstart(irrep1) + origStart1 + origIndex1,
                                                                    iHandler->getOrigNOCCstart(irrep2) + origStart2 + origIndex2,
                                                                    iHandler->getOrigNOCCstart(irrep3) + origStart3 + origIndex3,
                                                                    iHandler->getOrigNOCCstart(irrep4) + origStart4 + origIndex4 );
                                          }
                                       }
                                    }
                                 }
                                 
                                 //Rotate mem1_ajkl = U_ai mem2_ijkl (but only orig and target block1 indices)
                                 {
                                    char notra = 'N';
                                    double alpha = 1.0;
                                    double beta = 0.0; //SET !!!
                                    int rightdim = origSize2 * origSize3 * origSize4;
                                    int leftdim = targetSize1;
                                    int middledim = origSize1;
                                    double * rotationBlock = unitary->getBlock(irrep1) + targetStart1 + linsize1 * origStart1; // --> lda = linsize1;
                                    int lda = linsize1;
                                    dgemm_(&notra,&notra,&leftdim,&rightdim,&middledim, &alpha,rotationBlock,&lda,mem2,&middledim, &beta,mem1,&leftdim);
                                 }
                                 
                                 //Loop the target block 2
                                 int startTB2 = (irrep1==irrep2) ? targetBlock1 : 0;
                                 int stopTB2  = factorOA2 + factorV2;
                                 if (cutCorners){
                                    //case A: First index was forced to OA. Second one can be both OA and V
                                    //case B: First index can be OA and V. If the first one is V, the second one must be V as well.
                                    if ((caseB) && (!block1_OA)){ startTB2 = max( startTB2 , factorOA2 ); }
                                    //case C: First index can be OA and V. If the first index is V, the second one must be OA.
                                    if ((caseC) && (!block1_OA)){ stopTB2 = factorOA2; }
                                    //case D: At most two indices can be virtual: no restrictions yet.
                                 }
                                 for (int targetBlock2 = startTB2; targetBlock2 < stopTB2; targetBlock2++){
                     
                                    const bool block2_OA   = (targetBlock2 < factorOA2) ? true : false;
                                    const int targetStart2 = (block2_OA) ? (targetBlock2 * blocksizeOA2)
                                                                         : (linsizeOA2 + (targetBlock2 - factorOA2) * blocksizeV2);
                                    const int targetStop2  = (block2_OA) ? min( targetStart2 + blocksizeOA2 , linsizeOA2 )
                                                                         : min( targetStart2 + blocksizeV2  , linsize2   );
                                    const int targetSize2  = max( targetStop2 - targetStart2 , 0 );
                                    
                                    //Rotate mem2_abkl = U_bj mem1_ajkl (only block2 indices)
                                    {
                                       char trans = 'T';
                                       char notra = 'N';
                                       double alpha = 1.0;
                                       double beta = 0.0; //SET !!!
                                       int loop = origSize3 * origSize4;
                                       int jump_mem1 = targetSize1*origSize2;
                                       int jump_mem2 = targetSize1*targetSize2;
                                       int rightdim = targetSize2;
                                       int leftdim = targetSize1;
                                       int middledim = origSize2;
                                       double * rotationBlock = unitary->getBlock(irrep2) + targetStart2 + linsize2 * origStart2;  // --> lda = linsize2;
                                       int ldb = linsize2;
                                       for (int cntloop=0; cntloop<loop; cntloop++){
                                          dgemm_(&notra,&trans,&leftdim,&rightdim,&middledim, &alpha,mem1+cntloop*jump_mem1,&leftdim,rotationBlock,&ldb,
                                                                                              &beta, mem2+cntloop*jump_mem2,&leftdim);
                                       }
                                    }
                                    
                                    //Loop the target block 3
                                    int startTB3 = (irrep1==irrep3) ? targetBlock1 : 0;
                                    int stopTB3  = factorOA3 + factorV3;
                                    if (cutCorners){
                                       //All cases: at most two indices can be V. If the first two are V, the third must be OA.
                                       if ((!block1_OA) && (!block2_OA)){ stopTB3 = factorOA3; }
                                       //case A: First index is OA. Second index can be OA and V. If second index is V, the third one must be OA.
                                       if ((caseA) && (!block2_OA)){ stopTB3 = factorOA3; }
                                       //case B: First index V  --> second one V as well. Third must be OA, but imposed by first if statement above.
                                       //        First index OA --> (2,3) can be (OA,OA); (OA,V); (V,OA) and (V,V) --> no restrictions.
                                       if (caseC){
                                          if (!block1_OA){ startTB3 = max( startTB3 , factorOA3 ); } // First index V  --> (2,3) must be (OA,V)
                                          if ((block1_OA) && (!block2_OA)){ stopTB3 = factorOA3; } // First index OA --> if index2 is V, the index3 must be OA.
                                       }
                                    }
                                    for (int targetBlock3 = startTB3; targetBlock3 < stopTB3; targetBlock3++){
                        
                                       const bool block3_OA   = (targetBlock3 < factorOA3) ? true : false;
                                       const int targetStart3 = (block3_OA) ? (targetBlock3 * blocksizeOA3)
                                                                            : (linsizeOA3 + (targetBlock3 - factorOA3) * blocksizeV3);
                                       const int targetStop3  = (block3_OA) ? min( targetStart3 + blocksizeOA3 , linsizeOA3 )
                                                                            : min( targetStart3 + blocksizeV3  , linsize3   );
                                       const int targetSize3  = max( targetStop3 - targetStart3 , 0 );
                                       
                                       //Rotate mem3_abcl = U_ck mem2_abkl (only block3 indices)
                                       {
                                          char trans = 'T';
                                          char notra = 'N';
                                          double alpha = 1.0;
                                          double beta = 0.0; //SET !!!
                                          int loop = origSize4;
                                          int jump_mem2 = targetSize1*targetSize2*origSize3;
                                          int jump_mem3 = targetSize1*targetSize2*targetSize3;
                                          int rightdim = targetSize3;
                                          int leftdim = targetSize1*targetSize2;
                                          int middledim = origSize3;
                                          double * rotationBlock = unitary->getBlock(irrep3) + targetStart3 + linsize3 * origStart3;  // --> lda = linsize3;
                                          int ldb = linsize3;
                                          for (int cntloop=0; cntloop<loop; cntloop++){
                                             dgemm_(&notra,&trans,&leftdim,&rightdim,&middledim, &alpha,mem2+cntloop*jump_mem2,&leftdim,rotationBlock,&ldb,
                                                                                                 &beta, mem3+cntloop*jump_mem3,&leftdim);
                                          }
                                       }
                                       
                                       //Calculate V_partial_abcd = U_dl mem3_abcl and add to the FourIndexObject for all relevant d-indices
                                       const int leftdim = targetSize1 * targetSize2 * targetSize3;
                                       
                                       #pragma omp parallel for schedule(static)
                                       for (int counter=0; counter<leftdim; counter++){
                                          
                                          int targetIndex1 = counter % targetSize1;
                                          int temp = (counter - targetIndex1) / targetSize1;
                                          int targetIndex2 = temp % targetSize2;
                                          int targetIndex3 = (temp - targetIndex2) / targetSize2;
                                          
                                          const int hamIndex1 = iHandler->getOrigNOCCstart(irrep1) + targetStart1 + targetIndex1;
                                          const int hamIndex2 = iHandler->getOrigNOCCstart(irrep2) + targetStart2 + targetIndex2;
                                          const int hamIndex3 = iHandler->getOrigNOCCstart(irrep3) + targetStart3 + targetIndex3;
                                          
                                          //Only once per unique matrix element, does the relevant term need to be added
                                          if ((hamIndex1 <= hamIndex2) && (hamIndex1 <= hamIndex3)){
                                          
                                             //Loop the target INDEX 4
                                             int startTI4 = (irrep2==irrep4) ? targetStart2 + targetIndex2 : 0; //hamIndex2 <= hamIndex4
                                             int stopTI4  = linsizeOA4 + linsizeV4;
                                             if (cutCorners){
                                                //Get the number of virtuals
                                                int numV = 0;
                                                if (!block1_OA){ numV++; }
                                                if (!block2_OA){ numV++; }
                                                if (!block3_OA){ numV++; }
                                                //All cases: max two virtuals
                                                if (numV>=2){ stopTI4 = linsizeOA4; }
                                                //case A: With (1,2,3)=(OA,OA,OA) index4 can be both OA and V
                                                //        With (1,2,3)=(OA,OA,V) index4 van be both OA and V
                                                //        With (1,2,3)=(OA,V,OA) index4 must be V!
                                                if ((caseA) && (!block2_OA)){ startTI4 = max( startTI4, linsizeOA4 ); }
                                                //case B: For all cases not constrained by numV<=2, index4 can be both OA and V.
                                                //case C: For all cases not constrained by numV<=2, index4 can be both OA and V, except if the second index is V
                                                if ((caseC) && (!block2_OA)){ startTI4 = max( startTI4, linsizeOA4 ); }
                                             }
                                             for (int targetIndex4 = startTI4; targetIndex4 < stopTI4; targetIndex4++){
                                             
                                                const int hamIndex4  = iHandler->getOrigNOCCstart(irrep4) + targetIndex4;
                                             
                                                //Only once per unique matrix element, add contribution; hamIndex2 <= hamIndex4 ensured by targetStart4
                                                if ( (hamIndex1 != hamIndex2) || ( (hamIndex1 == hamIndex2) && (hamIndex3 <= hamIndex4) ) ){
                                                
                                                   double value = 0.0;
                                                   double * rotatedBlock = unitary->getBlock(irrep4) + linsize4 * origStart4;
                                                   for (int origIndex4=0; origIndex4<origSize4; origIndex4++){
                                                      value += mem3[counter + leftdim * origIndex4] * rotatedBlock[targetIndex4 + linsize4 * origIndex4];
                                                   }
                                                   VmatRotated->add(irrep1, irrep2, irrep3, irrep4, targetStart1 + targetIndex1,
                                                                    targetStart2 + targetIndex2, targetStart3 + targetIndex3, targetIndex4, value);
                                                   
                                                }
                                                
                                             }
                                          }
                                       }//End of of last rotation V_partial_abcd = U_dl mem3_abcl
                                    }
                                 }
                              }
                           } // End of loop original blocks
                        }
                     }
                  }
               } // End of non-zero irrep block (all linsizes > 0)
            }
         }
      }
   }

}

void CheMPS2::DMRGSCFVmatRotations::fillVmatDMRGBlockWise(Hamiltonian * HamDMRG, DMRGSCFunitary * unitary, double * mem1, double * mem2, double * mem3, const int maxBlockSize) const{
   
   /***********************************************************
   **   Two-body terms; use eightfold permutation symmetry   **
   **   Requires 3 buffers of size maxBlockSize^4            **
   **   maxBlockSize = ceil(maxlinsize/factor) at alloc      **
   ***********************************************************/
   for (int irrep1 = 0; irrep1<numberOfIrreps; irrep1++){
      for (int irrep2 = irrep1; irrep2<numberOfIrreps; irrep2++){
         const int productSymm = Irreps::directProd(irrep1,irrep2);
         for (int irrep3 = irrep1; irrep3<numberOfIrreps; irrep3++){
            const int irrep4 = Irreps::directProd(productSymm,irrep3);
            if (irrep4>=irrep2){
            
               const int linsizeDMRG1 = iHandler->getNDMRG(irrep1);
               const int linsizeDMRG2 = iHandler->getNDMRG(irrep2);
               const int linsizeDMRG3 = iHandler->getNDMRG(irrep3);
               const int linsizeDMRG4 = iHandler->getNDMRG(irrep4);
               
               int factorDMRG1 = max( (int) ( ceil((1.0 * linsizeDMRG1) / maxBlockSize) + 0.01 ) , 1 ); //factor >= linsize/maxBlockSize
               int factorDMRG2 = max( (int) ( ceil((1.0 * linsizeDMRG2) / maxBlockSize) + 0.01 ) , 1 );
               int factorDMRG3 = max( (int) ( ceil((1.0 * linsizeDMRG3) / maxBlockSize) + 0.01 ) , 1 );
               
               const int blocksizeDMRG1 = min( (int) ( ceil( (1.0 * linsizeDMRG1) / factorDMRG1 ) + 0.01 ) , maxBlockSize ); //Hence at most maxBlockSize
               const int blocksizeDMRG2 = min( (int) ( ceil( (1.0 * linsizeDMRG2) / factorDMRG2 ) + 0.01 ) , maxBlockSize );
               const int blocksizeDMRG3 = min( (int) ( ceil( (1.0 * linsizeDMRG3) / factorDMRG3 ) + 0.01 ) , maxBlockSize );
               
               while (factorDMRG1 * blocksizeDMRG1 < linsizeDMRG1){ factorDMRG1++; }
               while (factorDMRG2 * blocksizeDMRG2 < linsizeDMRG2){ factorDMRG2++; }
               while (factorDMRG3 * blocksizeDMRG3 < linsizeDMRG3){ factorDMRG3++; }
               
               if ((linsizeDMRG1>0) && (linsizeDMRG2>0) && (linsizeDMRG3>0) && (linsizeDMRG4>0)){
               
                  const int linsizeORIG1 = iHandler->getNORB(irrep1);
                  const int linsizeORIG2 = iHandler->getNORB(irrep2);
                  const int linsizeORIG3 = iHandler->getNORB(irrep3);
                  const int linsizeORIG4 = iHandler->getNORB(irrep4);
               
                  int factorORIG1 = max( (int) ( ceil((1.0 * linsizeORIG1) / maxBlockSize) + 0.01 ) , 1 ); //factor >= linsize/maxBlockSize
                  int factorORIG2 = max( (int) ( ceil((1.0 * linsizeORIG2) / maxBlockSize) + 0.01 ) , 1 );
                  int factorORIG3 = max( (int) ( ceil((1.0 * linsizeORIG3) / maxBlockSize) + 0.01 ) , 1 );
                  int factorORIG4 = max( (int) ( ceil((1.0 * linsizeORIG4) / maxBlockSize) + 0.01 ) , 1 );
                  
                  const int blocksizeORIG1 = min( (int) ( ceil( (1.0 * linsizeORIG1) / factorORIG1 ) + 0.01 ) , maxBlockSize ); //Hence at most maxBlockSize
                  const int blocksizeORIG2 = min( (int) ( ceil( (1.0 * linsizeORIG2) / factorORIG2 ) + 0.01 ) , maxBlockSize );
                  const int blocksizeORIG3 = min( (int) ( ceil( (1.0 * linsizeORIG3) / factorORIG3 ) + 0.01 ) , maxBlockSize );
                  const int blocksizeORIG4 = min( (int) ( ceil( (1.0 * linsizeORIG4) / factorORIG4 ) + 0.01 ) , maxBlockSize );
                  
                  while (factorORIG1 * blocksizeORIG1 < linsizeORIG1){ factorORIG1++; }
                  while (factorORIG2 * blocksizeORIG2 < linsizeORIG2){ factorORIG2++; }
                  while (factorORIG3 * blocksizeORIG3 < linsizeORIG3){ factorORIG3++; }
                  while (factorORIG4 * blocksizeORIG4 < linsizeORIG4){ factorORIG4++; }
                  
                  //Reset the FourIndex object
                  for (int cnt1=0; cnt1<linsizeDMRG1; cnt1++){
                     for (int cnt2=0; cnt2<linsizeDMRG2; cnt2++){
                        for (int cnt3=0; cnt3<linsizeDMRG3; cnt3++){
                           for (int cnt4=0; cnt4<linsizeDMRG4; cnt4++){
                              HamDMRG->setVmat( iHandler->getDMRGcumulative(irrep1) + cnt1, iHandler->getDMRGcumulative(irrep2) + cnt2,
                                                iHandler->getDMRGcumulative(irrep3) + cnt3, iHandler->getDMRGcumulative(irrep4) + cnt4, 0.0 );
                           }
                        }
                     }
                  }
                  
                  //Loop original blocks
                  for (int origBlock1=0; origBlock1<factorORIG1; origBlock1++){
                  
                     const int origStart1 = origBlock1 * blocksizeORIG1;
                     const int origStop1  = min( (origBlock1 + 1) * blocksizeORIG1 , linsizeORIG1 );
                     const int origSize1  = max( origStop1 - origStart1 , 0 );
                     
                     for (int origBlock2=0; origBlock2<factorORIG2; origBlock2++){
                     
                        const int origStart2 = origBlock2 * blocksizeORIG2;
                        const int origStop2  = min( (origBlock2 + 1) * blocksizeORIG2 , linsizeORIG2 );
                        const int origSize2  = max( origStop2 - origStart2 , 0 );
                     
                        for (int origBlock3=0; origBlock3<factorORIG3; origBlock3++){
                        
                           const int origStart3 = origBlock3 * blocksizeORIG3;
                           const int origStop3  = min( (origBlock3 + 1) * blocksizeORIG3 , linsizeORIG3 );
                           const int origSize3  = max( origStop3 - origStart3 , 0 );
                        
                           for (int origBlock4=0; origBlock4<factorORIG4; origBlock4++){
                           
                              const int origStart4 = origBlock4 * blocksizeORIG4;
                              const int origStop4  = min( (origBlock4 + 1) * blocksizeORIG4 , linsizeORIG4 );
                              const int origSize4  = max( origStop4 - origStart4 , 0 );
                              
                              //Loop the target block 1
                              for (int dmrgBlock1=0; dmrgBlock1<factorDMRG1; dmrgBlock1++){
                  
                                 const int dmrgStart1 = dmrgBlock1 * blocksizeDMRG1;
                                 const int dmrgStop1  = min( (dmrgBlock1 + 1) * blocksizeDMRG1 , linsizeDMRG1 );
                                 const int dmrgSize1  = max( dmrgStop1 - dmrgStart1 , 0 );
                                 
                                 //Copy HamOrig->getVmat(i,j,k,l) for the particular ORIGINAL block into mem2_ijkl
                                 for (int origIndex1=0; origIndex1<origSize1; origIndex1++){
                                    for (int origIndex2=0; origIndex2<origSize2; origIndex2++){
                                       for (int origIndex3=0; origIndex3<origSize3; origIndex3++){
                                          for (int origIndex4=0; origIndex4<origSize4; origIndex4++){
                                             mem2[ origIndex1 + origSize1 * (origIndex2 + origSize2 * (origIndex3 + origSize3 * origIndex4) ) ]
                                                = HamOrig->getVmat( iHandler->getOrigNOCCstart(irrep1) + origStart1 + origIndex1,
                                                                    iHandler->getOrigNOCCstart(irrep2) + origStart2 + origIndex2,
                                                                    iHandler->getOrigNOCCstart(irrep3) + origStart3 + origIndex3,
                                                                    iHandler->getOrigNOCCstart(irrep4) + origStart4 + origIndex4 );
                                          }
                                       }
                                    }
                                 }
                                 
                                 //Rotate mem1_ajkl = U_ai mem2_ijkl (but only orig and target block1 indices)
                                 {
                                    char notra = 'N';
                                    double alpha = 1.0;
                                    double beta  = 0.0; //SET !!!
                                    int rightdim  = origSize2 * origSize3 * origSize4;
                                    int leftdim   = dmrgSize1;
                                    int middledim = origSize1;
                                    double * rotationBlock = unitary->getBlock(irrep1) + iHandler->getNOCC(irrep1) + dmrgStart1 + linsizeORIG1 * origStart1;
                                    int lda = linsizeORIG1;
                                    dgemm_(&notra,&notra,&leftdim,&rightdim,&middledim, &alpha,rotationBlock,&lda,mem2,&middledim, &beta,mem1,&leftdim);
                                 }
                                 
                                 //Loop the target block 2
                                 for (int dmrgBlock2=((irrep1==irrep2)?dmrgBlock1:0); dmrgBlock2<factorDMRG2; dmrgBlock2++){
                     
                                    const int dmrgStart2 = dmrgBlock2 * blocksizeDMRG2;
                                    const int dmrgStop2  = min( (dmrgBlock2 + 1) * blocksizeDMRG2 , linsizeDMRG2 );
                                    const int dmrgSize2  = max( dmrgStop2 - dmrgStart2 , 0 );
                                    
                                    //Rotate mem2_abkl = U_bj mem1_ajkl (only block2 indices)
                                    {
                                       char trans = 'T';
                                       char notra = 'N';
                                       double alpha = 1.0;
                                       double beta  = 0.0; //SET !!!
                                       int loop = origSize3 * origSize4;
                                       int jump_mem1 = dmrgSize1 * origSize2;
                                       int jump_mem2 = dmrgSize1 * dmrgSize2;
                                       int rightdim  = dmrgSize2;
                                       int leftdim   = dmrgSize1;
                                       int middledim = origSize2;
                                       double * rotationBlock = unitary->getBlock(irrep2) + iHandler->getNOCC(irrep2) + dmrgStart2 + linsizeORIG2 * origStart2;
                                       int ldb = linsizeORIG2;
                                       for (int cntloop=0; cntloop<loop; cntloop++){
                                          dgemm_(&notra,&trans,&leftdim,&rightdim,&middledim, &alpha,mem1+cntloop*jump_mem1,&leftdim,rotationBlock,&ldb,
                                                                                              &beta, mem2+cntloop*jump_mem2,&leftdim);
                                       }
                                    }
                                    
                                    //Loop the target block 3
                                    for (int dmrgBlock3=((irrep1==irrep3)?dmrgBlock1:0); dmrgBlock3<factorDMRG3; dmrgBlock3++){
                        
                                       const int dmrgStart3 = dmrgBlock3 * blocksizeDMRG3;
                                       const int dmrgStop3 = min( (dmrgBlock3 + 1) * blocksizeDMRG3 , linsizeDMRG3 );
                                       const int dmrgSize3 = max( dmrgStop3 - dmrgStart3 , 0 );
                                       
                                       //Rotate mem3_abcl = U_ck mem2_abkl (only block3 indices)
                                       {
                                          char trans = 'T';
                                          char notra = 'N';
                                          double alpha = 1.0;
                                          double beta = 0.0; //SET !!!
                                          //int loop = origSize4;
                                          int jump_mem2 = dmrgSize1*dmrgSize2*origSize3;
                                          int jump_mem3 = dmrgSize1*dmrgSize2*dmrgSize3;
                                          int rightdim  = dmrgSize3;
                                          int leftdim   = dmrgSize1*dmrgSize2;
                                          int middledim = origSize3;
                                          double * rotationBlock = unitary->getBlock(irrep3) + iHandler->getNOCC(irrep3) + dmrgStart3 + linsizeORIG3 * origStart3;
                                          int ldb = linsizeORIG3;
                                          for (int cntloop=0; cntloop<origSize4; cntloop++){
                                             dgemm_(&notra,&trans,&leftdim,&rightdim,&middledim, &alpha,mem2+cntloop*jump_mem2,&leftdim,rotationBlock,&ldb,
                                                                                                 &beta, mem3+cntloop*jump_mem3,&leftdim);
                                          }
                                       }
                                       
                                       //Calculate V_partial_abcd = U_dl mem3_abcl and add to the FourIndexObject for all relevant d-indices
                                       const int leftdim = dmrgSize1 * dmrgSize2 * dmrgSize3;
                                       
                                       #pragma omp parallel for schedule(static)
                                       for (int counter=0; counter<leftdim; counter++){
                                          
                                          int dmrgIndex1 = counter % dmrgSize1;
                                          int temp       = (counter - dmrgIndex1) / dmrgSize1;
                                          int dmrgIndex2 = temp    % dmrgSize2;
                                          int dmrgIndex3 = (temp    - dmrgIndex2) / dmrgSize2;
                                          
                                          const int hamIndex1 = iHandler->getDMRGcumulative(irrep1) + dmrgStart1 + dmrgIndex1;
                                          const int hamIndex2 = iHandler->getDMRGcumulative(irrep2) + dmrgStart2 + dmrgIndex2;
                                          const int hamIndex3 = iHandler->getDMRGcumulative(irrep3) + dmrgStart3 + dmrgIndex3;
                                          
                                          //Only once per unique matrix element, does the relevant term need to be added
                                          if ((hamIndex1 <= hamIndex2) && (hamIndex1 <= hamIndex3)){
                                          
                                             int dmrgStart4 = (irrep2==irrep4) ? dmrgStart2 + dmrgIndex2 : 0; //hamIndex2 <= hamIndex4
                                          
                                             for (int dmrgIndex4=dmrgStart4; dmrgIndex4<linsizeDMRG4; dmrgIndex4++){
                                             
                                                const int hamIndex4 = iHandler->getDMRGcumulative(irrep4) + dmrgIndex4;
                                             
                                                //Only once per unique matrix element, add contribution; hamIndex2 <= hamIndex4 ensured by targetStart4
                                                if ( (hamIndex1 != hamIndex2) || ( (hamIndex1 == hamIndex2) && (hamIndex3 <= hamIndex4) ) ){
                                                
                                                   double value = 0.0;
                                                   double * rotatedBlock = unitary->getBlock(irrep4) + iHandler->getNOCC(irrep4) + linsizeORIG4 * origStart4;
                                                   for (int origIndex4=0; origIndex4<origSize4; origIndex4++){
                                                      value += mem3[counter + leftdim * origIndex4] * rotatedBlock[dmrgIndex4 + linsizeORIG4 * origIndex4];
                                                   }
                                                   HamDMRG->addToVmat( hamIndex1, hamIndex2, hamIndex3, hamIndex4, value );
                                                   
                                                }
                                                
                                             }
                                          }
                                       }//End of of last rotation V_partial_abcd = U_dl mem3_abcl
                                    }
                                 }
                              }
                           } // End of loop original blocks
                        }
                     }
                  }
               } // End of non-zero irrep block (all linsizes > 0)
            }
         }
      }
   }

}

void CheMPS2::DMRGSCFVmatRotations::fillRotatedTEIBlockWise(DMRGSCFintegrals * theRotatedTEI, DMRGSCFunitary * unitary, double * mem1, double * mem2, double * mem3, const int maxBlockSize) const{

   // First do Coulomb object : ( c1 <= c2 | a3 <= a4 )
   for (int Ic1 = 0; Ic1 < numberOfIrreps; Ic1++){
      for (int Ic2 = Ic1; Ic2 < numberOfIrreps; Ic2++){
         const int Icc = Irreps::directProd( Ic1, Ic2 );
         for (int Ia3 = 0; Ia3 < numberOfIrreps; Ia3++){
            const int Ia4 = Irreps::directProd( Ia3, Icc );
            if ( Ia3 <= Ia4 ){
               
               int linsize1 = iHandler->getNORB( Ic1 );
               int linsize2 = iHandler->getNORB( Ic2 );
               int linsize3 = iHandler->getNORB( Ia3 );
               int linsize4 = iHandler->getNORB( Ia4 );
               
               int linsizeOA1 = iHandler->getNOCC( Ic1 ) + iHandler->getNDMRG( Ic1 );
               int linsizeOA2 = iHandler->getNOCC( Ic2 ) + iHandler->getNDMRG( Ic2 );
               
               if (( linsizeOA1 > 0 ) && ( linsizeOA2 > 0 ) && ( linsize3 > 0 ) && ( linsize4 > 0 )){
               
                  int factor1 = max( (int) ( ceil((1.0 * linsize1) / maxBlockSize) + 0.01 ) , 1 ); //factor >= linsize/maxBlockSize
                  int factor2 = max( (int) ( ceil((1.0 * linsize2) / maxBlockSize) + 0.01 ) , 1 );
                  int factor3 = max( (int) ( ceil((1.0 * linsize3) / maxBlockSize) + 0.01 ) , 1 );
                  int factor4 = max( (int) ( ceil((1.0 * linsize4) / maxBlockSize) + 0.01 ) , 1 );
                  
                  const int blocksize1 = min( (int) ( ceil( (1.0 * linsize1) / factor1 ) + 0.01 ) , maxBlockSize ); //Hence at most maxBlockSize
                  const int blocksize2 = min( (int) ( ceil( (1.0 * linsize2) / factor2 ) + 0.01 ) , maxBlockSize );
                  const int blocksize3 = min( (int) ( ceil( (1.0 * linsize3) / factor3 ) + 0.01 ) , maxBlockSize );
                  const int blocksize4 = min( (int) ( ceil( (1.0 * linsize4) / factor4 ) + 0.01 ) , maxBlockSize );
                  
                  while (factor1 * blocksize1 < linsize1){ factor1++; }
                  while (factor2 * blocksize2 < linsize2){ factor2++; }
                  while (factor3 * blocksize3 < linsize3){ factor3++; }
                  while (factor4 * blocksize4 < linsize4){ factor4++; }
                  
                  int factorOA1 = (linsizeOA1 == 0) ? 0 : max( (int) ( ceil((1.0 * linsizeOA1) / maxBlockSize) + 0.01 ) , 1 );
                  int factorOA2 = (linsizeOA2 == 0) ? 0 : max( (int) ( ceil((1.0 * linsizeOA2) / maxBlockSize) + 0.01 ) , 1 );
                  
                  const int blocksizeOA1 = (linsizeOA1 == 0) ? 1 : min( (int) ( ceil( (1.0 * linsizeOA1) / factorOA1 ) + 0.01 ) , maxBlockSize );
                  const int blocksizeOA2 = (linsizeOA2 == 0) ? 1 : min( (int) ( ceil( (1.0 * linsizeOA2) / factorOA2 ) + 0.01 ) , maxBlockSize );
                  
                  if (linsizeOA1 > 0){ while (factorOA1 * blocksizeOA1 < linsizeOA1){ factorOA1++; } }
                  if (linsizeOA2 > 0){ while (factorOA2 * blocksizeOA2 < linsizeOA2){ factorOA2++; } }
                  
                  //Clear the Coulomb object
                  for (int c1 = 0; c1 < linsizeOA1; c1++){
                     for (int c2 = 0; c2 < linsizeOA2; c2++){
                        for (int a3 = 0; a3 < linsize3; a3++){
                           for (int a4 = 0; a4 < linsize4; a4++){
                              theRotatedTEI->set_coulomb( Ic1, Ic2, Ia3, Ia4, c1, c2, a3, a4, 0.0 );
                           }
                        }
                     }
                  }
                  
                  //Loop original blocks
                  for (int origBlock1 = 0; origBlock1 < factor1; origBlock1++){
                  
                     const int origStart1 = origBlock1 * blocksize1;
                     const int origStop1  = min( (origBlock1 + 1) * blocksize1 , linsize1 );
                     const int origSize1  = max( origStop1 - origStart1 , 0 );
                     
                     for (int origBlock2 = 0; origBlock2 < factor2; origBlock2++){
                     
                        const int origStart2 = origBlock2 * blocksize2;
                        const int origStop2  = min( (origBlock2 + 1) * blocksize2 , linsize2 );
                        const int origSize2  = max( origStop2 - origStart2 , 0 );
                     
                        for (int origBlock3 = 0; origBlock3 < factor3; origBlock3++){
                        
                           const int origStart3 = origBlock3 * blocksize3;
                           const int origStop3  = min( (origBlock3 + 1) * blocksize3 , linsize3 );
                           const int origSize3  = max( origStop3 - origStart3 , 0 );
                        
                           for (int origBlock4 = 0; origBlock4 < factor4; origBlock4++){
                           
                              const int origStart4 = origBlock4 * blocksize4;
                              const int origStop4  = min( (origBlock4 + 1) * blocksize4 , linsize4 );
                              const int origSize4  = max( origStop4 - origStart4 , 0 );
                              
                              //Loop target blocks
                              for (int targetBlock1 = 0; targetBlock1 < factorOA1; targetBlock1++){
                              
                                 const int targetStart1 = targetBlock1 * blocksizeOA1;
                                 const int targetStop1  = min( (targetBlock1 + 1) * blocksizeOA1 , linsizeOA1 );
                                 const int targetSize1  = max( targetStop1 - targetStart1 , 0 );
                              
                                 //Copy HamOrig->getVmat(c1,a3,c2,a4) for the particular ORIGINAL block into mem2[c1,c2,a3,a4]
                                 for (int origIndex1 = 0; origIndex1 < origSize1; origIndex1++){
                                    for (int origIndex2 = 0; origIndex2 < origSize2; origIndex2++){
                                       for (int origIndex3 = 0; origIndex3 < origSize3; origIndex3++){
                                          for (int origIndex4 = 0; origIndex4 < origSize4; origIndex4++){
                                             mem2[ origIndex1 + origSize1 * (origIndex2 + origSize2 * (origIndex3 + origSize3 * origIndex4) ) ]
                                                = HamOrig->getVmat( iHandler->getOrigNOCCstart( Ic1 ) + origStart1 + origIndex1,
                                                                    iHandler->getOrigNOCCstart( Ia3 ) + origStart3 + origIndex3,
                                                                    iHandler->getOrigNOCCstart( Ic2 ) + origStart2 + origIndex2,
                                                                    iHandler->getOrigNOCCstart( Ia4 ) + origStart4 + origIndex4 );
                                          }
                                       }
                                    }
                                 }
                                 
                                 //Rotate mem1[ c1, j, k, l ] = U[ c1, i ] * mem2[ i, j, k, l ]
                                 {
                                    char notrans = 'N';
                                    double alpha = 1.0;
                                    double beta = 0.0; //SET !!!
                                    int rightdim = origSize2 * origSize3 * origSize4;
                                    int leftdim = targetSize1;
                                    int middledim = origSize1;
                                    double * rotationBlock = unitary->getBlock( Ic1 ) + targetStart1 + linsize1 * origStart1; // --> lda = linsize1;
                                    int lda = linsize1;
                                    dgemm_(&notrans, &notrans ,&leftdim, &rightdim, &middledim, &alpha,
                                           rotationBlock, &lda, mem2, &middledim, &beta, mem1, &leftdim);
                                 }
                                 
                                 //Loop target block 2
                                 for (int targetBlock2 = (( Icc==0 ) ? targetBlock1 : 0); targetBlock2 < factorOA2; targetBlock2++){
                                 
                                    const int targetStart2 = targetBlock2 * blocksizeOA2;
                                    const int targetStop2  = min( (targetBlock2 + 1) * blocksizeOA2 , linsizeOA2 );
                                    const int targetSize2  = max( targetStop2 - targetStart2 , 0 );
                                    
                                    //Rotate mem2[ c1, c2, k, l ] = U[ c2, j ] * mem1[ c1, j, k, l ]
                                    {
                                       char trans = 'T';
                                       char notrans = 'N';
                                       double alpha = 1.0;
                                       double beta = 0.0; //SET !!!
                                       int loopsize = origSize3 * origSize4;
                                       int jump_mem1 = targetSize1 * origSize2;
                                       int jump_mem2 = targetSize1 * targetSize2;
                                       int rightdim = targetSize2;
                                       int leftdim = targetSize1;
                                       int middledim = origSize2;
                                       double * rotationBlock = unitary->getBlock( Ic2 ) + targetStart2 + linsize2 * origStart2; // --> ldb = linsize2;
                                       int ldb = linsize2;
                                       for (int cntloop = 0; cntloop < loopsize; cntloop++){
                                          dgemm_(&notrans, &trans, &leftdim, &rightdim, &middledim, &alpha,
                                                 mem1 + cntloop * jump_mem1, &leftdim, rotationBlock, &ldb, &beta, mem2 + cntloop * jump_mem2, &leftdim);
                                       }
                                    }
                                    
                                    //Loop target block 3
                                    for (int targetBlock3 = 0; targetBlock3 < factor3; targetBlock3++){
                                    
                                       const int targetStart3 = targetBlock3 * blocksize3;
                                       const int targetStop3  = min( (targetBlock3 + 1) * blocksize3 , linsize3 );
                                       const int targetSize3  = max( targetStop3 - targetStart3 , 0 );
                                       
                                       //Rotate mem3[ c1, c2, a3, l ] = U[ a3, k ] * mem2[ c1, c2, l, l ]
                                       {
                                          char trans = 'T';
                                          char notrans = 'N';
                                          double alpha = 1.0;
                                          double beta = 0.0; //SET !!!
                                          int jump_mem2 = targetSize1 * targetSize2 * origSize3;
                                          int jump_mem3 = targetSize1 * targetSize2 * targetSize3;
                                          int rightdim = targetSize3;
                                          int leftdim = targetSize1 * targetSize2;
                                          int middledim = origSize3;
                                          double * rotationBlock = unitary->getBlock( Ia3 ) + targetStart3 + linsize3 * origStart3; // --> ldb = linsize3;
                                          int ldb = linsize3;
                                          for (int cntloop = 0; cntloop < origSize4; cntloop++){
                                             dgemm_(&notrans, &trans, &leftdim, &rightdim, &middledim, &alpha,
                                                    mem2 + cntloop * jump_mem2, &leftdim, rotationBlock, &ldb, &beta, mem3 + cntloop * jump_mem3, &leftdim);
                                          }
                                       }
                                       
                                       //Calculate ( c1 c2 | a3 a4 )_partial and add to the relevant parts of the Coulomb object
                                       const int loopsize = targetSize1 * targetSize2 * targetSize3;
                                       
                                       #pragma omp parallel for schedule(static)
                                       for (int counter = 0; counter < loopsize; counter++){
                                       
                                          const int c1_rel = counter % targetSize1;
                                          int temp = ( counter - c1_rel ) / targetSize1;
                                          const int c2_rel = temp % targetSize2;
                                          const int a3_rel = ( temp - c2_rel ) / targetSize2;
                                          
                                          const int c1 = c1_rel + targetStart1;
                                          const int c2 = c2_rel + targetStart2;
                                          const int a3 = a3_rel + targetStart3;
                                          
                                          int a4start      = 0;
                                          const int a4stop = linsize4;
                                          // If Icc==0 --> be careful that ( c1 <= c2 | a3 <= a4 ) valid
                                          // If Icc!=0 --> automatically OK because Ic1 < Ic2 and Ia3 < Ia4 checked at beginning
                                          if ( Icc == 0 ){ a4start = (( c1 <= c2 ) ? a3 : a4stop); }
                                          for (int a4 = a4start; a4 < a4stop; a4++){
                                             
                                             double value = 0.0;
                                             double * rotatedBlock = unitary->getBlock( Ia4 ) + linsize4 * origStart4;
                                             for (int origIndex4 = 0; origIndex4 < origSize4; origIndex4++){
                                                value += mem3[counter + loopsize * origIndex4] * rotatedBlock[a4 + linsize4 * origIndex4];
                                             }
                                             
                                             theRotatedTEI->add_coulomb( Ic1, Ic2, Ia3, Ia4, c1, c2, a3, a4, value );
                                             
                                          }
                                       }
                                    }
                                 }
                              }//targetBlock1
                           }
                        }
                     }
                  }//origBlock1
               }
            }
         }
      }
   }
   
   // Then do Exchange object : ( c1 v3 | c2 v4 )
   for (int Ic1 = 0; Ic1 < numberOfIrreps; Ic1++){
      for (int Ic2 = Ic1; Ic2 < numberOfIrreps; Ic2++){
         const int Icc = Irreps::directProd( Ic1, Ic2 );
         for (int Iv3 = 0; Iv3 < numberOfIrreps; Iv3++){
            const int Iv4 = Irreps::directProd( Iv3, Icc ); //No restriction on Iv4, only Ic1 <= Ic2
            
            int linsize1 = iHandler->getNORB( Ic1 );
            int linsize2 = iHandler->getNORB( Ic2 );
            int linsize3 = iHandler->getNORB( Iv3 );
            int linsize4 = iHandler->getNORB( Iv4 );
            
            int linsizeOA1    = iHandler->getNOCC(  Ic1 ) + iHandler->getNDMRG( Ic1 );
            int linsizeOA2    = iHandler->getNOCC(  Ic2 ) + iHandler->getNDMRG( Ic2 );
            int linsizeV3     = iHandler->getNVIRT( Iv3 );
            int linsizeV4     = iHandler->getNVIRT( Iv4 );
            const int shiftv3 = iHandler->getNOCC(  Iv3 ) + iHandler->getNDMRG( Iv3 );
            const int shiftv4 = iHandler->getNOCC(  Iv4 ) + iHandler->getNDMRG( Iv4 );
            
            if (( linsizeOA1 > 0 ) && ( linsizeOA2 > 0 ) && ( linsizeV3 > 0 ) && ( linsizeV4 > 0 )){
            
               int factor1 = max( (int) ( ceil((1.0 * linsize1) / maxBlockSize) + 0.01 ) , 1 ); //factor >= linsize/maxBlockSize
               int factor2 = max( (int) ( ceil((1.0 * linsize2) / maxBlockSize) + 0.01 ) , 1 );
               int factor3 = max( (int) ( ceil((1.0 * linsize3) / maxBlockSize) + 0.01 ) , 1 );
               int factor4 = max( (int) ( ceil((1.0 * linsize4) / maxBlockSize) + 0.01 ) , 1 );
               
               const int blocksize1 = min( (int) ( ceil( (1.0 * linsize1) / factor1 ) + 0.01 ) , maxBlockSize ); //Hence at most maxBlockSize
               const int blocksize2 = min( (int) ( ceil( (1.0 * linsize2) / factor2 ) + 0.01 ) , maxBlockSize );
               const int blocksize3 = min( (int) ( ceil( (1.0 * linsize3) / factor3 ) + 0.01 ) , maxBlockSize );
               const int blocksize4 = min( (int) ( ceil( (1.0 * linsize4) / factor4 ) + 0.01 ) , maxBlockSize );
               
               while (factor1 * blocksize1 < linsize1){ factor1++; }
               while (factor2 * blocksize2 < linsize2){ factor2++; }
               while (factor3 * blocksize3 < linsize3){ factor3++; }
               while (factor4 * blocksize4 < linsize4){ factor4++; }
               
               int factorOA1 = (linsizeOA1 == 0) ? 0 : max( (int) ( ceil((1.0 * linsizeOA1) / maxBlockSize) + 0.01 ) , 1 );
               int factorOA2 = (linsizeOA2 == 0) ? 0 : max( (int) ( ceil((1.0 * linsizeOA2) / maxBlockSize) + 0.01 ) , 1 );
               int factorV3  = (linsizeV3  == 0) ? 0 : max( (int) ( ceil((1.0 * linsizeV3 ) / maxBlockSize) + 0.01 ) , 1 );
               int factorV4  = (linsizeV4  == 0) ? 0 : max( (int) ( ceil((1.0 * linsizeV4 ) / maxBlockSize) + 0.01 ) , 1 );
               
               const int blocksizeOA1 = (linsizeOA1 == 0) ? 1 : min( (int) ( ceil( (1.0 * linsizeOA1) / factorOA1 ) + 0.01 ) , maxBlockSize );
               const int blocksizeOA2 = (linsizeOA2 == 0) ? 1 : min( (int) ( ceil( (1.0 * linsizeOA2) / factorOA2 ) + 0.01 ) , maxBlockSize );
               const int blocksizeV3  = (linsizeV3  == 0) ? 1 : min( (int) ( ceil( (1.0 * linsizeV3 ) / factorV3  ) + 0.01 ) , maxBlockSize );
               const int blocksizeV4  = (linsizeV4  == 0) ? 1 : min( (int) ( ceil( (1.0 * linsizeV4 ) / factorV4  ) + 0.01 ) , maxBlockSize );
               
               if (linsizeOA1 > 0){ while (factorOA1 * blocksizeOA1 < linsizeOA1){ factorOA1++; } }
               if (linsizeOA2 > 0){ while (factorOA2 * blocksizeOA2 < linsizeOA2){ factorOA2++; } }
               if (linsizeV3  > 0){ while (factorV3  * blocksizeV3  < linsizeV3 ){ factorV3++;  } }
               if (linsizeV4  > 0){ while (factorV4  * blocksizeV4  < linsizeV4 ){ factorV4++;  } }
               
               //Clear the Exchange object
               for (int c1 = 0; c1 < linsizeOA1; c1++){
                  for (int c2 = 0; c2 < linsizeOA2; c2++){
                     for (int v3 = 0; v3 < linsizeV3; v3++){
                        for (int v4 = 0; v4 < linsizeV4; v4++){
                           theRotatedTEI->set_exchange( Ic1, Ic2, Iv3, Iv4, c1, c2, shiftv3 + v3, shiftv4 + v4, 0.0 );
                        }
                     }
                  }
               }
               
               //Loop original blocks
               for (int origBlock1 = 0; origBlock1 < factor1; origBlock1++){
               
                  const int origStart1 = origBlock1 * blocksize1;
                  const int origStop1  = min( (origBlock1 + 1) * blocksize1 , linsize1 );
                  const int origSize1  = max( origStop1 - origStart1 , 0 );
                  
                  for (int origBlock2 = 0; origBlock2 < factor2; origBlock2++){
                  
                     const int origStart2 = origBlock2 * blocksize2;
                     const int origStop2  = min( (origBlock2 + 1) * blocksize2 , linsize2 );
                     const int origSize2  = max( origStop2 - origStart2 , 0 );
                  
                     for (int origBlock3 = 0; origBlock3 < factor3; origBlock3++){
                     
                        const int origStart3 = origBlock3 * blocksize3;
                        const int origStop3  = min( (origBlock3 + 1) * blocksize3 , linsize3 );
                        const int origSize3  = max( origStop3 - origStart3 , 0 );
                     
                        for (int origBlock4 = 0; origBlock4 < factor4; origBlock4++){
                        
                           const int origStart4 = origBlock4 * blocksize4;
                           const int origStop4  = min( (origBlock4 + 1) * blocksize4 , linsize4 );
                           const int origSize4  = max( origStop4 - origStart4 , 0 );
                           
                           //Loop target blocks
                           for (int targetBlock1 = 0; targetBlock1 < factorOA1; targetBlock1++){
                           
                              const int targetStart1 = targetBlock1 * blocksizeOA1;
                              const int targetStop1  = min( (targetBlock1 + 1) * blocksizeOA1 , linsizeOA1 );
                              const int targetSize1  = max( targetStop1 - targetStart1 , 0 );
                           
                              //Copy HamOrig->getVmat(c1,c2,v3,v4) for the particular ORIGINAL block into mem2[c1,c2,v3,v4]
                              for (int origIndex1 = 0; origIndex1 < origSize1; origIndex1++){
                                 for (int origIndex2 = 0; origIndex2 < origSize2; origIndex2++){
                                    for (int origIndex3 = 0; origIndex3 < origSize3; origIndex3++){
                                       for (int origIndex4 = 0; origIndex4 < origSize4; origIndex4++){
                                          mem2[ origIndex1 + origSize1 * (origIndex2 + origSize2 * (origIndex3 + origSize3 * origIndex4) ) ]
                                             = HamOrig->getVmat( iHandler->getOrigNOCCstart( Ic1 ) + origStart1 + origIndex1,
                                                                 iHandler->getOrigNOCCstart( Ic2 ) + origStart2 + origIndex2,
                                                                 iHandler->getOrigNOCCstart( Iv3 ) + origStart3 + origIndex3,
                                                                 iHandler->getOrigNOCCstart( Iv4 ) + origStart4 + origIndex4 );
                                       }
                                    }
                                 }
                              }
                              
                              //Rotate mem1[ c1, j, k, l ] = U[ c1, i ] * mem2[ i, j, k, l ]
                              {
                                 char notrans = 'N';
                                 double alpha = 1.0;
                                 double beta = 0.0; //SET !!!
                                 int rightdim = origSize2 * origSize3 * origSize4;
                                 int leftdim = targetSize1;
                                 int middledim = origSize1;
                                 double * rotationBlock = unitary->getBlock( Ic1 ) + targetStart1 + linsize1 * origStart1; // --> lda = linsize1;
                                 int lda = linsize1;
                                 dgemm_(&notrans, &notrans ,&leftdim, &rightdim, &middledim, &alpha,
                                        rotationBlock, &lda, mem2, &middledim, &beta, mem1, &leftdim);
                              }
                              
                              //Loop target block 2
                              for (int targetBlock2 = (( Icc==0 ) ? targetBlock1 : 0); targetBlock2 < factorOA2; targetBlock2++){
                              
                                 const int targetStart2 = targetBlock2 * blocksizeOA2;
                                 const int targetStop2  = min( (targetBlock2 + 1) * blocksizeOA2 , linsizeOA2 );
                                 const int targetSize2  = max( targetStop2 - targetStart2 , 0 );
                                 
                                 //Rotate mem2[ c1, c2, k, l ] = U[ c2, j ] * mem1[ c1, j, k, l ]
                                 {
                                    char trans = 'T';
                                    char notrans = 'N';
                                    double alpha = 1.0;
                                    double beta = 0.0; //SET !!!
                                    int loopsize = origSize3 * origSize4;
                                    int jump_mem1 = targetSize1 * origSize2;
                                    int jump_mem2 = targetSize1 * targetSize2;
                                    int rightdim = targetSize2;
                                    int leftdim = targetSize1;
                                    int middledim = origSize2;
                                    double * rotationBlock = unitary->getBlock( Ic2 ) + targetStart2 + linsize2 * origStart2; // --> ldb = linsize2;
                                    int ldb = linsize2;
                                    for (int cntloop = 0; cntloop < loopsize; cntloop++){
                                       dgemm_(&notrans, &trans, &leftdim, &rightdim, &middledim, &alpha,
                                              mem1 + cntloop * jump_mem1, &leftdim, rotationBlock, &ldb, &beta, mem2 + cntloop * jump_mem2, &leftdim);
                                    }
                                 }
                                 
                                 //Loop target block 3
                                 for (int targetBlock3 = 0; targetBlock3 < factorV3; targetBlock3++){
                                 
                                    const int targetStart3 = shiftv3 + targetBlock3 * blocksizeV3;
                                    const int targetStop3  = shiftv3 + min( (targetBlock3 + 1) * blocksizeV3 , linsizeV3 );
                                    const int targetSize3  = max( targetStop3 - targetStart3 , 0 );
                                    
                                    //Rotate mem3[ c1, c2, v3, l ] = U[ v3, k ] * mem2[ c1, c2, l, l ]
                                    {
                                       char trans = 'T';
                                       char notrans = 'N';
                                       double alpha = 1.0;
                                       double beta = 0.0; //SET !!!
                                       int jump_mem2 = targetSize1 * targetSize2 * origSize3;
                                       int jump_mem3 = targetSize1 * targetSize2 * targetSize3;
                                       int rightdim = targetSize3;
                                       int leftdim = targetSize1 * targetSize2;
                                       int middledim = origSize3;
                                       double * rotationBlock = unitary->getBlock( Iv3 ) + targetStart3 + linsize3 * origStart3; // --> ldb = linsize3;
                                       int ldb = linsize3;
                                       for (int cntloop = 0; cntloop < origSize4; cntloop++){
                                          dgemm_(&notrans, &trans, &leftdim, &rightdim, &middledim, &alpha,
                                                 mem2 + cntloop * jump_mem2, &leftdim, rotationBlock, &ldb, &beta, mem3 + cntloop * jump_mem3, &leftdim);
                                       }
                                    }
                                    
                                    //Calculate ( c1 c2 | v3 v4 )_partial and add to the relevant parts of the Exchange object
                                    const int loopsize = targetSize1 * targetSize2 * targetSize3;
                                    
                                    #pragma omp parallel for schedule(static)
                                    for (int counter = 0; counter < loopsize; counter++){
                                    
                                       const int c1_rel = counter % targetSize1;
                                       int temp = ( counter - c1_rel ) / targetSize1;
                                       const int c2_rel = temp % targetSize2;
                                       const int v3_rel = ( temp - c2_rel ) / targetSize2;
                                       
                                       const int c1 = c1_rel + targetStart1;
                                       const int c2 = c2_rel + targetStart2;
                                       const int v3 = v3_rel + targetStart3;
                                       
                                       int v4start      = shiftv4;
                                       const int v4stop = linsize4;
                                       // If Icc==0 --> be careful that c1 <= c2 is valid
                                       // If Icc!=0 --> automatically OK because Ic1 < Ic2
                                       if (( Icc == 0 ) && ( c1 > c2 )){ v4start = v4stop; }
                                       for (int v4 = v4start; v4 < v4stop; v4++){
                                          
                                          double value = 0.0;
                                          double * rotatedBlock = unitary->getBlock( Iv4 ) + linsize4 * origStart4;
                                          for (int origIndex4 = 0; origIndex4 < origSize4; origIndex4++){
                                             value += mem3[counter + loopsize * origIndex4] * rotatedBlock[v4 + linsize4 * origIndex4];
                                          }
                                          
                                          theRotatedTEI->add_exchange( Ic1, Ic2, Iv3, Iv4, c1, c2, v3, v4, value );
                                          
                                       }
                                    }
                                 }
                              }
                           }//targetBlock1
                        }
                     }
                  }
               }//origBlock1
            }
         }
      }
   }

}


