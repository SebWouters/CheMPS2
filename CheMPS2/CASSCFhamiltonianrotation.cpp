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

#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include "CASSCF.h"
#include "Lapack.h"

using std::min;
using std::max;

void CheMPS2::CASSCF::fillRotatedHamAllInMemory(double * temp1, double * temp2){

   //Constant part of the energy
   HamRotated->setEconst(HamOrig->getEconst());
   
   //One-body terms: diagonal in the irreps.
   int passed = 0;
   for (int irrep=0; irrep<numberOfIrreps; irrep++){
   
      int linsize = iHandler->getNORB(irrep);
      if (linsize>1){
         
         for (int cnt1=0; cnt1<linsize; cnt1++){
            for (int cnt2=0; cnt2<linsize; cnt2++){
               temp1[cnt1 + linsize * cnt2] = HamOrig->getTmat(passed+cnt1,passed+cnt2);
            }
         }
         
         char trans = 'T';
         char notra = 'N';
         double alpha = 1.0;
         double beta = 0.0;
         dgemm_(&notra,&notra,&linsize,&linsize,&linsize,&alpha, unitary->getBlock(irrep),&linsize,temp1,&linsize,&beta,temp2,&linsize);
         dgemm_(&notra,&trans,&linsize,&linsize,&linsize,&alpha,temp2,&linsize, unitary->getBlock(irrep),&linsize,&beta,temp1,&linsize);
         
         for (int cnt1=0; cnt1<linsize; cnt1++){
            for (int cnt2=cnt1; cnt2<linsize; cnt2++){
               HamRotated->setTmat(passed+cnt1,passed+cnt2, temp1[cnt1 + linsize * cnt2] ) ;
            }
         }
         
      }
      if (linsize==1){ HamRotated->setTmat(passed,passed, HamOrig->getTmat(passed,passed) ); }
      
      passed += linsize;
   
   }
   
   //Two-body terms --> use eightfold permutation symmetry
   for (int irrep1 = 0; irrep1<numberOfIrreps; irrep1++){
      for (int irrep2 = irrep1; irrep2<numberOfIrreps; irrep2++){
         const int productSymm = SymmInfo.directProd(irrep1,irrep2);
         for (int irrep3 = irrep1; irrep3<numberOfIrreps; irrep3++){
            const int irrep4 = SymmInfo.directProd(productSymm,irrep3);
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
                              HamRotated->setVmat( iHandler->getOrigNOCCstart(irrep1) + cnt1, iHandler->getOrigNOCCstart(irrep2) + cnt2,
                                                   iHandler->getOrigNOCCstart(irrep3) + cnt3, iHandler->getOrigNOCCstart(irrep4) + cnt4,
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

void CheMPS2::CASSCF::fillRotatedHamInMemoryBlockWise(double * mem1, double * mem2, double * mem3, const int maxBlockSize){

   /************************************
   **   Constant part of the energy   **
   ************************************/
   HamRotated->setEconst(HamOrig->getEconst());
   
   /**********************************************************
   **   One-body terms; diagonal in the irreps.             **
   **   Requires mem1 and mem2 to be of size maxlinsize^2   **
   **********************************************************/
   int passed = 0;
   for (int irrep=0; irrep<numberOfIrreps; irrep++){
   
      int linsize = iHandler->getNORB(irrep);
      if (linsize>1){
         
         for (int cnt1=0; cnt1<linsize; cnt1++){
            for (int cnt2=0; cnt2<linsize; cnt2++){
               mem1[cnt1 + linsize * cnt2] = HamOrig->getTmat(passed+cnt1,passed+cnt2);
            }
         }
         
         char trans = 'T';
         char notra = 'N';
         double alpha = 1.0;
         double beta = 0.0;
         dgemm_(&notra,&notra,&linsize,&linsize,&linsize,&alpha, unitary->getBlock(irrep),&linsize,mem1,&linsize,&beta,mem2,&linsize);
         dgemm_(&notra,&trans,&linsize,&linsize,&linsize,&alpha,mem2,&linsize, unitary->getBlock(irrep),&linsize,&beta,mem1,&linsize);
         
         for (int cnt1=0; cnt1<linsize; cnt1++){
            for (int cnt2=cnt1; cnt2<linsize; cnt2++){
               HamRotated->setTmat(passed+cnt1,passed+cnt2, mem1[cnt1 + linsize * cnt2] ) ;
            }
         }
         
      }
      if (linsize==1){ HamRotated->setTmat(passed,passed, HamOrig->getTmat(passed,passed) ); }
      passed += linsize;
   
   }
   
   /***********************************************************
   **   Two-body terms; use eightfold permutation symmetry   **
   **   Requires 3 buffers of size maxBlockSize^4            **
   **   maxBlockSize = ceil(maxlinsize/factor) at alloc      **
   ***********************************************************/
   for (int irrep1 = 0; irrep1<numberOfIrreps; irrep1++){
      for (int irrep2 = irrep1; irrep2<numberOfIrreps; irrep2++){
         const int productSymm = SymmInfo.directProd(irrep1,irrep2);
         for (int irrep3 = irrep1; irrep3<numberOfIrreps; irrep3++){
            const int irrep4 = SymmInfo.directProd(productSymm,irrep3);
            if (irrep4>=irrep2){
            
               const int linsize1 = iHandler->getNORB(irrep1);
               const int linsize2 = iHandler->getNORB(irrep2);
               const int linsize3 = iHandler->getNORB(irrep3);
               const int linsize4 = iHandler->getNORB(irrep4);
               
               if ((linsize1>0) && (linsize2>0) && (linsize3>0) && (linsize4>0)){
               
                  /** The scheme: 
                      Split up each of the linsizes into factor parts of size (linsize/factor). Make sure the boundaries are clear.
                      
                      Algorithm:
                         Reset the FourIndex object for the particular symmetry case under consideration --> HamRotated->setVmat(i,j,k,l,0.0)
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
                  
                  //Split up the linsizes
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
                  
                  /* for (int blocki=0; blocki<factori; blocki++)
                        int starti = blocki * blocksizei
                        int stopi = min( (blocki+1) * blocksizei , linsizei )
                        int sizei = upperboundi - starti
                        for (int indexi = starti; indexi<stopi; indexi++) 
                  */
                  
                  //Reset the FourIndex object
                  for (int cnt1=0; cnt1<linsize1; cnt1++){
                     for (int cnt2=0; cnt2<linsize2; cnt2++){
                        for (int cnt3=0; cnt3<linsize3; cnt3++){
                           for (int cnt4=0; cnt4<linsize4; cnt4++){
                              HamRotated->setVmat( iHandler->getOrigNOCCstart(irrep1) + cnt1, iHandler->getOrigNOCCstart(irrep2) + cnt2,
                                                   iHandler->getOrigNOCCstart(irrep3) + cnt3, iHandler->getOrigNOCCstart(irrep4) + cnt4, 0.0);
                           }
                        }
                     }
                  }
                  
                  //Loop original blocks
                  for (int origBlock1=0; origBlock1<factor1; origBlock1++){
                  
                     const int origStart1 = origBlock1 * blocksize1;
                     const int origStop1 = min( (origBlock1 + 1) * blocksize1 , linsize1 );
                     const int origSize1 = max( origStop1 - origStart1 , 0 );
                     
                     for (int origBlock2=0; origBlock2<factor2; origBlock2++){
                     
                        const int origStart2 = origBlock2 * blocksize2;
                        const int origStop2 = min( (origBlock2 + 1) * blocksize2 , linsize2 );
                        const int origSize2 = max( origStop2 - origStart2 , 0 );
                     
                        for (int origBlock3=0; origBlock3<factor3; origBlock3++){
                        
                           const int origStart3 = origBlock3 * blocksize3;
                           const int origStop3 = min( (origBlock3 + 1) * blocksize3 , linsize3 );
                           const int origSize3 = max( origStop3 - origStart3 , 0 );
                        
                           for (int origBlock4=0; origBlock4<factor4; origBlock4++){
                           
                              const int origStart4 = origBlock4 * blocksize4;
                              const int origStop4 = min( (origBlock4 + 1) * blocksize4 , linsize4 );
                              const int origSize4 = max( origStop4 - origStart4 , 0 );
                              
                              //Loop the target block 1
                              for (int targetBlock1=0; targetBlock1<factor1; targetBlock1++){
                  
                                 const int targetStart1 = targetBlock1 * blocksize1;
                                 const int targetStop1 = min( (targetBlock1 + 1) * blocksize1 , linsize1 );
                                 const int targetSize1 = max( targetStop1 - targetStart1 , 0 );
                                 
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
                                 for (int targetBlock2=((irrep1==irrep2)?targetBlock1:0); targetBlock2<factor2; targetBlock2++){
                     
                                    const int targetStart2 = targetBlock2 * blocksize2;
                                    const int targetStop2 = min( (targetBlock2 + 1) * blocksize2 , linsize2 );
                                    const int targetSize2 = max( targetStop2 - targetStart2 , 0 );
                                    
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
                                    for (int targetBlock3=((irrep1==irrep3)?targetBlock1:0); targetBlock3<factor3; targetBlock3++){
                        
                                       const int targetStart3 = targetBlock3 * blocksize3;
                                       const int targetStop3 = min( (targetBlock3 + 1) * blocksize3 , linsize3 );
                                       const int targetSize3 = max( targetStop3 - targetStart3 , 0 );
                                       
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
                                          
                                          const int rowindex  = targetIndex1 + targetSize1 * (targetIndex2 + targetSize2 * targetIndex3);
                                          const int hamIndex1 = iHandler->getOrigNOCCstart(irrep1) + targetStart1 + targetIndex1;
                                          const int hamIndex2 = iHandler->getOrigNOCCstart(irrep2) + targetStart2 + targetIndex2;
                                          const int hamIndex3 = iHandler->getOrigNOCCstart(irrep3) + targetStart3 + targetIndex3;
                                          
                                          //Only once per unique matrix element, does the relevant term need to be added
                                          if ((hamIndex1 <= hamIndex2) && (hamIndex1 <= hamIndex3)){
                                          
                                             int targetStart4 = (irrep2==irrep4) ? targetStart2 + targetIndex2 : 0; //hamIndex2 <= hamIndex4
                                          
                                             for (int targetIndex4=targetStart4; targetIndex4<linsize4; targetIndex4++){
                                             
                                                const int hamIndex4 = iHandler->getOrigNOCCstart(irrep4) + targetIndex4;
                                             
                                                //Only once per unique matrix element, add contribution; hamIndex2 <= hamIndex4 ensured by targetStart4
                                                if ( (hamIndex1 != hamIndex2) || ( (hamIndex1 == hamIndex2) && (hamIndex3 <= hamIndex4) ) ){
                                                
                                                   double value = 0.0;
                                                   double * rotatedBlock = unitary->getBlock(irrep4) + linsize4 * origStart4;
                                                   for (int origIndex4=0; origIndex4<origSize4; origIndex4++){
                                                      value += mem3[rowindex + leftdim * origIndex4] * rotatedBlock[targetIndex4 + linsize4 * origIndex4];
                                                   }
                                                   HamRotated->addToVmat( hamIndex1, hamIndex2, hamIndex3, hamIndex4, value );
                                                   
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


