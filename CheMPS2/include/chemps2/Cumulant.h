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

#ifndef CUMULANT_CHEMPS2_H
#define CUMULANT_CHEMPS2_H

#include "ThreeDM.h"
#include "TwoDM.h"

namespace CheMPS2{
/** Cumulant class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date November 26, 2015
    
    The cumulant class contains routines to approximate the spinfree 4-RDM \f$ \Gamma^4 \f$ by neglecting the fourth order cumulant \f$ \Lambda^4 \f$. Based on the spinfree density matrices
    \f{eqnarray*}{
      \Gamma^1_{ip}        & = & \sum\limits_{\sigma}                   \braket{ 0 | \hat{a}_{i \sigma}^{\dagger} \hat{a}_{p \sigma} | 0 } \\
      \Gamma^2_{ijpq}      & = & \sum\limits_{\sigma \tau}              \braket{ 0 | \hat{a}_{i \sigma}^{\dagger} \hat{a}_{j \tau}^{\dagger} \hat{a}_{q \tau} \hat{a}_{p \sigma} | 0 } \\
      \Gamma^3_{ijkpqr}    & = & \sum\limits_{\sigma \tau \chi}         \braket{ 0 | \hat{a}_{i \sigma}^{\dagger} \hat{a}_{j \tau}^{\dagger} \hat{a}_{k \chi}^{\dagger} \hat{a}_{r \chi} \hat{a}_{q \tau} \hat{a}_{p \sigma} | 0 } \\
      \Gamma^4_{ijklpqrs}  & = & \sum\limits_{\sigma \tau \chi \omega } \braket{ 0 | \hat{a}_{i \sigma}^{\dagger} \hat{a}_{j \tau}^{\dagger} \hat{a}_{k \chi}^{\dagger} \hat{a}_{l \omega}^{\dagger} \hat{a}_{s \omega} \hat{a}_{r \chi} \hat{a}_{q \tau} \hat{a}_{p \sigma} | 0 }
    \f}
    and the spinfree second order cumulant 
    \f{eqnarray*}{
     \Lambda^2_{ijpq} & = & \Gamma^2_{ijpq} - \Gamma^1_{ip} \Gamma^1_{jq} + \frac{1}{2} \Gamma^1_{iq} \Gamma^1_{jp}
    \f}
    the spinfree 4-RDM can be written as [CUM1]:
    \f{eqnarray*}{
                     &   &  \Gamma^4_{ijklpqrs} \\
                     & = &  \Lambda^4_{ijklpqrs} \\
                     & + &  \Gamma^3_{ijkpqr} \Gamma^1_{ls}
                                 - \frac{1}{2} \Gamma^3_{ijksqr} \Gamma^1_{lp}
                                 - \frac{1}{2} \Gamma^3_{ijkpsr} \Gamma^1_{lq}
                                 - \frac{1}{2} \Gamma^3_{ijkpqs} \Gamma^1_{lr} \\
                     & + &  \Gamma^3_{ijlpqs} \Gamma^1_{kr}
                                 - \frac{1}{2} \Gamma^3_{ijlrqs} \Gamma^1_{kp}
                                 - \frac{1}{2} \Gamma^3_{ijlprs} \Gamma^1_{kq}
                                 - \frac{1}{2} \Gamma^3_{ijlpqr} \Gamma^1_{ks} \\
                     & + &  \Gamma^3_{iklprs} \Gamma^1_{jq}
                                 - \frac{1}{2} \Gamma^3_{iklqrs} \Gamma^1_{jp}
                                 - \frac{1}{2} \Gamma^3_{iklpqs} \Gamma^1_{jr}
                                 - \frac{1}{2} \Gamma^3_{iklprq} \Gamma^1_{js} \\
                     & + &  \Gamma^3_{jklqrs} \Gamma^1_{ip}
                                 - \frac{1}{2} \Gamma^3_{jklprs} \Gamma^1_{iq}
                                 - \frac{1}{2} \Gamma^3_{jklqps} \Gamma^1_{ir}
                                 - \frac{1}{2} \Gamma^3_{jklqrp} \Gamma^1_{is} \\
                     & - &             \Gamma^2_{ijpq} \Gamma^2_{klrs}
                                 + \frac{1}{2} \Gamma^2_{ijpr} \Gamma^2_{klqs}
                                 + \frac{1}{2} \Gamma^2_{ijps} \Gamma^2_{klrq}
                                 + \frac{1}{2} \Gamma^2_{ijrq} \Gamma^2_{klps}
                                 + \frac{1}{2} \Gamma^2_{ijsq} \Gamma^2_{klrp} \\
                     & - & \frac{1}{3} \Gamma^2_{ijrs} \Gamma^2_{klpq}
                                 - \frac{1}{6} \Gamma^2_{ijrs} \Gamma^2_{klqp}
                                 - \frac{1}{6} \Gamma^2_{ijsr} \Gamma^2_{klpq}
                                 - \frac{1}{3} \Gamma^2_{ijsr} \Gamma^2_{klqp} \\
                     & - &             \Gamma^2_{ikpr} \Gamma^2_{jlqs}
                                 + \frac{1}{2} \Gamma^2_{ikpq} \Gamma^2_{jlrs}
                                 + \frac{1}{2} \Gamma^2_{ikps} \Gamma^2_{jlqr}
                                 + \frac{1}{2} \Gamma^2_{ikqr} \Gamma^2_{jlps}
                                 + \frac{1}{2} \Gamma^2_{iksr} \Gamma^2_{jlqp} \\
                     & - & \frac{1}{3} \Gamma^2_{ikqs} \Gamma^2_{jlpr}
                                 - \frac{1}{6} \Gamma^2_{iksq} \Gamma^2_{jlpr}
                                 - \frac{1}{6} \Gamma^2_{ikqs} \Gamma^2_{jlrp}
                                 - \frac{1}{3} \Gamma^2_{iksq} \Gamma^2_{jlrp} \\
                     & - &             \Gamma^2_{ilps} \Gamma^2_{kjrq}
                                 + \frac{1}{2} \Gamma^2_{ilpr} \Gamma^2_{kjsq}
                                 + \frac{1}{2} \Gamma^2_{ilpq} \Gamma^2_{kjrs}
                                 + \frac{1}{2} \Gamma^2_{ilrs} \Gamma^2_{kjpq}
                                 + \frac{1}{2} \Gamma^2_{ilqs} \Gamma^2_{kjrp} \\
                     & - & \frac{1}{3} \Gamma^2_{ilrq} \Gamma^2_{kjps}
                                 - \frac{1}{6} \Gamma^2_{ilqr} \Gamma^2_{kjps}
                                 - \frac{1}{6} \Gamma^2_{ilrq} \Gamma^2_{kjsp}
                                 - \frac{1}{3} \Gamma^2_{ilqr} \Gamma^2_{kjsp} \\
                     & + & 2           \Lambda^2_{ijpq} \Lambda^2_{klrs}
                                 -             \Lambda^2_{ijpr} \Lambda^2_{klqs}
                                 -             \Lambda^2_{ijps} \Lambda^2_{klrq}
                                 -             \Lambda^2_{ijrq} \Lambda^2_{klps}
                                 -             \Lambda^2_{ijsq} \Lambda^2_{klrp} \\
                     & + & \frac{2}{3} \Lambda^2_{ijrs} \Lambda^2_{klpq}
                                 + \frac{1}{3} \Lambda^2_{ijrs} \Lambda^2_{klqp}
                                 + \frac{1}{3} \Lambda^2_{ijsr} \Lambda^2_{klpq}
                                 + \frac{2}{3} \Lambda^2_{ijsr} \Lambda^2_{klqp} \\
                     & + & 2           \Lambda^2_{ikpr} \Lambda^2_{jlqs}
                                 -             \Lambda^2_{ikpq} \Lambda^2_{jlrs}
                                 -             \Lambda^2_{ikps} \Lambda^2_{jlqr}
                                 -             \Lambda^2_{ikqr} \Lambda^2_{jlps}
                                 -             \Lambda^2_{iksr} \Lambda^2_{jlqp} \\
                     & + & \frac{2}{3} \Lambda^2_{ikqs} \Lambda^2_{jlpr}
                                 + \frac{1}{3} \Lambda^2_{iksq} \Lambda^2_{jlpr}
                                 + \frac{1}{3} \Lambda^2_{ikqs} \Lambda^2_{jlrp}
                                 + \frac{2}{3} \Lambda^2_{iksq} \Lambda^2_{jlrp} \\
                     & + & 2           \Lambda^2_{ilps} \Lambda^2_{kjrq}
                                 -             \Lambda^2_{ilpr} \Lambda^2_{kjsq}
                                 -             \Lambda^2_{ilpq} \Lambda^2_{kjrs}
                                 -             \Lambda^2_{ilrs} \Lambda^2_{kjpq}
                                 -             \Lambda^2_{ilqs} \Lambda^2_{kjrp} \\
                     & + & \frac{2}{3} \Lambda^2_{ilrq} \Lambda^2_{kjps}
                                 + \frac{1}{3} \Lambda^2_{ilqr} \Lambda^2_{kjps}
                                 + \frac{1}{3} \Lambda^2_{ilrq} \Lambda^2_{kjsp}
                                 + \frac{2}{3} \Lambda^2_{ilqr} \Lambda^2_{kjsp}
    \f}
    By neglecting \f$ \Lambda^4 \f$, the cumulant approximation of the 4-RDM \f$ \Gamma^4 \f$ is obtained. \n
    \n
    [CUM1] M. Saitow, Y. Kurashige and T. Yanai, Journal of Chemical Physics 139, 044118 (2013). http://dx.doi.org/10.1063/1.4816627 \n*/
   class Cumulant{

      public:

         //! Get the cumulant approximation of \f$ \Gamma^4_{ijklpqrs} \f$, using HAM indices
         /** \param prob Pointer to the DMRG problem
             \param the3DM Pointer to the DMRG 3-RDM
             \param the2DM Pointer to the DMRG 2-RDM
             \param i index 1 of \f$ \Gamma^4_{ijklpqrs} \f$
             \param j index 2 of \f$ \Gamma^4_{ijklpqrs} \f$
             \param k index 3 of \f$ \Gamma^4_{ijklpqrs} \f$
             \param l index 4 of \f$ \Gamma^4_{ijklpqrs} \f$
             \param p index 5 of \f$ \Gamma^4_{ijklpqrs} \f$
             \param q index 6 of \f$ \Gamma^4_{ijklpqrs} \f$
             \param r index 7 of \f$ \Gamma^4_{ijklpqrs} \f$
             \param s index 8 of \f$ \Gamma^4_{ijklpqrs} \f$
             \return the desired value */
         static double gamma4_ham(const Problem * prob, const ThreeDM * the3DM, const TwoDM * the2DM, const int i, const int j, const int k, const int l, const int p, const int q, const int r, const int s);
         
         //! Contract the CASPT2 Fock operator with the cumulant approximation of \f$ \Gamma^4 \f$ in \f$ \mathcal{O}(L^7) \f$ time, using HAM indices
         /** \param prob Pointer to the DMRG problem
             \param the3DM Pointer to the DMRG 3-RDM
             \param the2DM Pointer to the DMRG 2-RDM
             \param fock Contains the SYMMETRIC fock operator \f$ F_{ls} \f$ = fock[l+L*s] = fock[s+L*l]
             \param result Contains the contraction: result[i+L*(j+L*(k+L*(p+L*(q+L*r))))] = \f$ \sum\limits_{ls} F_{ls} \Gamma^4_{ijklpqrs} \f$ */
         static void gamma4_fock_contract_ham(const Problem * prob, const ThreeDM * the3DM, const TwoDM * the2DM, double * fock, double * result);
         
      private:
      
         // Get the second order cumulant \f$ \Lambda^2_{ijpq} \f$, using HAM indices
         static double lambda2_ham(const TwoDM * the2DM, const int i, const int j, const int p, const int q);
         
   };
}

#endif
