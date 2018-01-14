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

#ifndef DMRGSCFWTILDE_CHEMPS2_H
#define DMRGSCFWTILDE_CHEMPS2_H

#include "DMRGSCFindices.h"

namespace CheMPS2{
/** DMRGSCF w_tilde class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date January 26, 2015
    
    Container class for the tensor \f$\tilde{w}_{pqrs}\f$. The definition and context can be found in CASSCF.h. For convenience, I repeat the definition of \f$\tilde{w}_{pqrs}\f$ here (remember that \f$I_p = I_q\f$ and \f$I_r = I_s\f$):
    \f{eqnarray*}{
      (p,r) & \in & (occ,occ) : \tilde{w}_{pqrs} \\
                          & = & 4 \delta_{pr}^{occ} \left[ h_{qs} + Q^{occ}_{qs} + Q^{act}_{qs} \right] \\
                          & + & 4 \left[ 4 (qp | sr) - ( qs | pr ) - ( qr | sp ) \right] \\
      (p,r) & \in & (act,act) : \tilde{w}_{pqrs} = 2 \Gamma^{1,act}_{rp} \left[ h_{qs} + Q^{occ}_{qs} \right] \\
                          & + & 2 \sum\limits_{\alpha\beta \in act} \left[ \Gamma^{2A,act}_{r \alpha p \beta} (qs | \alpha \beta ) + \left( \Gamma^{2A,act}_{r \alpha \beta p} + \Gamma^{2A,act}_{r p \beta \alpha} \right) (q \alpha | s \beta ) \right] \\
      (p,r) & \in & (act,occ) : \tilde{w}_{pqrs} \\
                          & = & 2 \sum\limits_{\alpha \in act} \Gamma^{1,act}_{\alpha p} \left[ 4 (q \alpha | s r) - (qs | \alpha r) - (qr | s \alpha) \right] \\
      (p,r) & \in & (occ,act) : \tilde{w}_{pqrs} \\
                          & = & 2 \sum\limits_{\beta \in act} \Gamma^{1,act}_{r \beta} \left[ 4 (q p | s \beta) -  (qs | p \beta) - (q \beta | sp) \right]
    \f}
*/
   class DMRGSCFwtilde{

      public:
      
         //! Constructor
         /** \param iHandler_in The DMRGSCFindices which contain information on the occupied, active, and virtual spaces */
         DMRGSCFwtilde(DMRGSCFindices * iHandler_in);
         
         //! Destructor
         virtual ~DMRGSCFwtilde();
         
         //! Clear
         void clear();
         
         //! Set an element of w_tilde_pqrs
         /** \param irrep_pq The irrep number of the first two indices pq
             \param irrep_rs The irrep number of the last two indices rs
             \param p The first index (within the symmetry block)
             \param q The second index (within the symmetry block)
             \param r The third index (within the symmetry block)
             \param s The fourth index (within the symmetry block)
             \param val The value to which the element of the tensor should be set */
         void set(const int irrep_pq, const int irrep_rs, const int p, const int q, const int r, const int s, const double val);

         //! Get an element of w_tilde_pqrs
         /** \param irrep_pq The irrep number of the first two indices pq
             \param irrep_rs The irrep number of the last two indices rs
             \param p The first index (within the symmetry block)
             \param q The second index (within the symmetry block)
             \param r The third index (within the symmetry block)
             \param s The fourth index (within the symmetry block)
             \return The requested element */
         double get(const int irrep_pq, const int irrep_rs, const int p, const int q, const int r, const int s) const;
         
         //! Get the (pr) subblock of w_tilde_pqrs, which is stored as w_tilde[ I_pq ][ I_rs ][ p + ( Nocc[I_pq] + Ndmrg[I_pq] ) * r ][ q + Ntotal[I_pq] * s ]
         /** \param irrep_pq The irrep number of the first two indices pq
             \param irrep_rs The irrep number of the last two indices rs
             \param p The first index (within the symmetry block)
             \param r The third index (within the symmetry block)
             \return Pointer to the requested subblock */
         double * getBlock(const int irrep_pq, const int irrep_rs, const int p, const int r);
      
      private:
      
         // The information on the occupied, active, and virtual spaces
         DMRGSCFindices * iHandler;
         
         int * Nocc_dmrg;
         
         // The elements: w_tilde[ I_pq ][ I_rs ][ p + ( Nocc[I_pq] + Ndmrg[I_pq] ) * r ][ q + Ntotal[I_pq] * s ]
         double **** wmattilde;

   };
}

#endif

