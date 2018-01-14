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

#ifndef IRREPS_CHEMPS2_H
#define IRREPS_CHEMPS2_H

#include <string>
using std::string;

namespace CheMPS2{
/** Irreps class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date February 7, 2013
    
    This class contains the symmetry group and irrep conventions. The program requires Abelian point groups with real character tables, with hence \f$I_{\alpha} \otimes I_{\alpha} = I_{trivial}\f$.
  
    \section irreps_conv Irrep conventions
    
    The same conventions as in Psi4 (beta5) are used. For convenience, they are listed below:\n
    \latexonly
    \begin{tabular}{|l|cccccccc|}
    \hline
      Symmetry Conventions & \multicolumn{8}{c|}{ Irrep Number \& Name } \\
    \hline
      Group Number \& Name & 0  & 1   & 2   & 3   & 4  & 5   & 6   & 7   \\
    \hline
      0: c1                & A  &     &     &     &    &     &     &     \\
      1: ci                & Ag & Au  &     &     &    &     &     &     \\
      2: c2                & A  & B   &     &     &    &     &     &     \\
      3: cs                & A' & A'' &     &     &    &     &     &     \\
      4: d2                & A  & B1  & B2  & B3  &    &     &     &     \\
      5: c2v               & A1 & A2  & B1  & B2  &    &     &     &     \\
      6: c2h               & Ag & Bg  & Au  & Bu  &    &     &     &     \\
      7: d2h               & Ag & B1g & B2g & B3g & Au & B1u & B2u & B3u \\
    \hline
    \end{tabular}
    \endlatexonly
    \htmlonly
    <table border="1">
    <tr><td> Symmetry Conventions </td><td colspan="8"> Irrep Number & Name </td></tr>
    <tr><td> Group Number & Name </td><td> 0  </td><td> 1   </td><td> 2   </td><td> 3   </td><td> 4  </td><td> 5   </td><td> 6   </td><td> 7   </td></tr>
    <tr><td> 0: c1               </td><td> A  </td><td>     </td><td>     </td><td>     </td><td>    </td><td>     </td><td>     </td><td>     </td></tr>
    <tr><td> 1: ci               </td><td> Ag </td><td> Au  </td><td>     </td><td>     </td><td>    </td><td>     </td><td>     </td><td>     </td></tr>
    <tr><td> 2: c2               </td><td> A  </td><td> B   </td><td>     </td><td>     </td><td>    </td><td>     </td><td>     </td><td>     </td></tr>
    <tr><td> 3: cs               </td><td> A' </td><td> A'' </td><td>     </td><td>     </td><td>    </td><td>     </td><td>     </td><td>     </td></tr>
    <tr><td> 4: d2               </td><td> A  </td><td> B1  </td><td> B2  </td><td> B3  </td><td>    </td><td>     </td><td>     </td><td>     </td></tr>
    <tr><td> 5: c2v              </td><td> A1 </td><td> A2  </td><td> B1  </td><td> B2  </td><td>    </td><td>     </td><td>     </td><td>     </td></tr>
    <tr><td> 6: c2h              </td><td> Ag </td><td> Bg  </td><td> Au  </td><td> Bu  </td><td>    </td><td>     </td><td>     </td><td>     </td></tr>
    <tr><td> 7: d2h              </td><td> Ag </td><td> B1g </td><td> B2g </td><td> B3g </td><td> Au </td><td> B1u </td><td> B2u </td><td> B3u </td></tr>
    </table>
    \endhtmlonly
    Note that these conventions allow to use the XOR operation for irrep multiplication. */
   class Irreps{

      public:
      
         //! Constructor
         Irreps();
         
         //! Constructor 2
         /** \param nGroup The group number (0 <= nGroup <= 7; else isActivated remains false) */
         Irreps(const int nGroup);
         
         //! Destructor
         virtual ~Irreps();
         
         //! Set the group
         /** \param nGroup Number from 0 to 7 (7 included)
             \return Validity of the group number (error returns false) */
         bool setGroup(const int nGroup);
         
         //! Whether the group number is already activated
         /** \return Whether the group number is already activated */
         bool getIsActivated() const;
         
         //! Get the group number
         /** \return The group number (-1 means not activated) */
         int getGroupNumber() const;
         
         //! Get the name of the group
         /** \return The group name ("error" means not activated) */
         string getGroupName() const;
         
         //! Get the name of the group corresponding to nGroup
         /** \param nGroup Group number
             \return The group name corresponding to nGroup */
         static string getGroupName(const int nGroup);
         
         //! Get the number of irreps for the currently activated group
         /** \return The number of irreps for the currently activated group (-1 means not activated) */
         int getNumberOfIrreps() const;
         
         //! Get the number of irreps for a certain group number
         /** \param nGroup The group number for which the number of irreps will be returned
             \return The number of irreps of that group (-1 means wrong group number) */
         static int getNumberOfIrreps(const int nGroup);
         
         //! Get the name of the irrep with number irrepNumber of the activated group. The irrep with number 0 is always the trivial irrep.
         /** \param irrepNumber The irrep number
             \return The irrep name (not activated returns "error1"; wrong number returns "error2") */
         string getIrrepName(const int irrepNumber) const;
         
         //! Get the direct product of the irreps with numbers Irrep1 and Irrep2: a bitwise XOR for psi4's conventions
         /** \param Irrep1 The number of the first irrep
             \param Irrep2 The number of the second irrep
             \return The direct product of I1 and I2 */
         static int directProd(const int Irrep1, const int Irrep2){ return Irrep1 ^ Irrep2; }
         
         //! Fill the array psi2molpro with the irrep conventions of molpro for the currently activated group
         /** \param psi2molpro The array to be filled: psi2molpro[psi_irrep] = molpro_irrep */
         void symm_psi2molpro( int * psi2molpro ) const;
         
         //! Fill the array psi2molpro with the irrep conventions of molpro for the group with symmetry label SymmLabel
         /** \param psi2molpro The array to be filled: psi2molpro[psi_irrep] = molpro_irrep
             \param SymmLabel The group for which psi2molpro needs to be filled */
         static void symm_psi2molpro( int * psi2molpro, const string SymmLabel );
         
         //! Print all info contained in this class
         static void printAll();
      
      private:
      
         //whether a relevant group number is currently set
         bool isActivated;
         
         //the currently set group number (check isActivated)
         int groupNumber;
         
         //number of irreps
         int nIrreps;
         
         //static member functions containing the group names, the number of irreps, and the names of the irreps
         static string getGroupNamePrivate(const int nGroup);
         static string getIrrepNamePrivate(const int nGroup, const int nIrrep);

   };
}

#endif

