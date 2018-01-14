.. CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
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

.. index:: DMRG algorithm

DMRG algorithm
==============

The density matrix renormalization group (DMRG) was first used for ab initio quantum chemistry in 1999 [DMRG1]_. The method variationally optimizes a low-rank tensor approximation of the full configuration interaction (FCI) solution. Suppose we have :math:`L` spatial orbitals. The FCI solution can in general be written as

.. math::

    \left|\Psi\right\rangle & = & \sum\limits_{\{ n_{i\sigma} \}} C^{n_{1\uparrow} n_{1\downarrow} n_{2\uparrow} n_{2\downarrow} n_{3\uparrow} ... n_{L\uparrow} n_{L\downarrow} }  \left( \hat{a}^{\dagger}_{1\uparrow} \right)^{n_{1\uparrow}} \left( \hat{a}^{\dagger}_{1\downarrow} \right)^{n_{1\downarrow}} \left( \hat{a}^{\dagger}_{2\uparrow} \right)^{n_{2\uparrow}} ... \left( \hat{a}^{\dagger}_{L\uparrow} \right)^{n_{L\uparrow}} \left( \hat{a}^{\dagger}_{L\downarrow} \right)^{n_{L\downarrow}} \left|-\right\rangle \nonumber \\
    & = & \sum\limits_{\{ n_{i\sigma} \}} C^{n_{1\uparrow} n_{1\downarrow} n_{2\uparrow} n_{2\downarrow} n_{3\uparrow} ... n_{L\uparrow} n_{L\downarrow} } \left| n_{1\uparrow} n_{1\downarrow} n_{2\uparrow} n_{2\downarrow} n_{3\uparrow} ... n_{L\uparrow} n_{L\downarrow} \right\rangle.

With successive singular value decompositions, the FCI :math:`C`-tensor can be composed into a matrix product state (MPS):

.. math::

    C^{ n_{1\uparrow} n_{1\downarrow} n_{2\uparrow} n_{2\downarrow} n_{3\uparrow} ... n_{L\uparrow} n_{L\downarrow} } & = & \sum\limits_{\alpha_1, \alpha_2, ..., \alpha_L} A[1]^{ n_{1\uparrow} n_{1\downarrow} }_{ \alpha_1 } A[2]^{n_{2\uparrow} n_{2\downarrow}}_{ \alpha_1 ; \alpha_2 } A[3]^{n_{3\uparrow} n_{3\downarrow}}_{ \alpha_2 ; \alpha_3 } ... A[L-1]^{n_{L-1\uparrow} n_{L-1\downarrow}}_{ \alpha_{L-2} ; \alpha_{L-1} } A[L]^{n_{L\uparrow} n_{L\downarrow}}_{ \alpha_{L-1} } \nonumber \\
    & = & \mathbf{A}[1]^{ n_{1\uparrow} n_{1\downarrow} } \mathbf{A}[2]^{n_{2\uparrow} n_{2\downarrow}} \mathbf{A}[3]^{n_{3\uparrow} n_{3\downarrow}} ... \mathbf{A}[L-1]^{n_{L-1\uparrow} n_{L-1\downarrow}} \mathbf{A}[L]^{n_{L\uparrow} n_{L\downarrow}},

where :math:`dim(\alpha_i) = min(4^i,4^{L-i})`. To make the method of polynomial complexity, the rank of the decomposition is truncated to a fixed dimension :math:`D`:

.. math::

   dim(\alpha_i) = min(4^i,4^{L-i},D).
   
The integer :math:`D` is called the bond, virtual, or auxiliary dimension. The DMRG algorithm consists of consecutive sweeps over the chain of orbitals, during which two neighbouring MPS tensors are variatonally optimized. Thereto they are combined into a two-orbital tensor:

.. math::

    \mathbf{B}[i]^{n_{i\uparrow} n_{i\downarrow} n_{i+1\uparrow} n_{i+1\downarrow}} = \mathbf{A}[i]^{n_{i\uparrow} n_{i\downarrow}} \mathbf{A}[i+1]^{n_{i+1\uparrow} n_{i+1\downarrow}}.

The Lagrangian

.. math::

    \mathcal{L} = \left\langle \Psi( \mathbf{B}[i] ) \mid \hat{H} \mid \Psi( \mathbf{B}[i] ) \right\rangle - E \left\langle \Psi( \mathbf{B}[i] ) \mid \Psi( \mathbf{B}[i] ) \right\rangle
    
is varied with respect to :math:`\mathbf{B}[i]` to yield an effective Hamiltonian eigenvalue equation. By exploiting the gauge freedom in the MPS, this eigenvalue equation can always be turned into a numerically stable standard eigenvalue equation for each local optimization step:

.. math::

    \mathbf{H}^{\text{effective}}[i] \times \mathbf{B}[i] = E \mathbf{B}[i].
    
Once :math:`\mathbf{B}[i]` is found, it is decomposed with a singular value decomposition:

.. math::

    B[i]^{n_{i\uparrow} n_{i\downarrow} n_{i+1\uparrow} n_{i+1\downarrow}}_{\alpha;\beta} = M_{(\alpha n_{i\uparrow} n_{i\downarrow});(n_{i+1\uparrow} n_{i+1\downarrow} \beta)} = \sum\limits_{\kappa} U_{(\alpha n_{i\uparrow} n_{i\downarrow});\kappa} \lambda_{\kappa} V^{\dagger}_{\kappa;(n_{i+1\uparrow} n_{i+1\downarrow} \beta)} = \sum\limits_{\kappa} A[i]^{n_{i\uparrow} n_{i\downarrow}}_{\alpha;\kappa} \lambda_{\kappa} A[i+1]^{n_{i+1\uparrow} n_{i+1\downarrow}}_{\kappa;\beta}.
    
For a normalized wavefunction :math:`\sum\limits_{\kappa} \lambda_{\kappa}^2 = 1`. If :math:`dim(\alpha) = dim(\beta) = D` then :math:`dim(\kappa) = 4D`. In order to keep the virtual dimension fixed to :math:`D`, the summation over :math:`\kappa` is truncated to the :math:`D` largest values :math:`\lambda_{\kappa}`. The discarded weight :math:`w_D[i] = \sum\limits_{\kappa > D} \lambda_{\kappa}^2` is a measure for the information loss.

For more information on the DMRG method, please read Ref. [DMRG2]_.

.. [DMRG1] S.R. White and R.L. Martin, *Journal of Chemical Physics* **110**, 4127 (1999), doi: `10.1063/1.478295 <http://dx.doi.org/10.1063/1.478295>`_
.. [DMRG2] S. Wouters and D. Van Neck, *European Physical Journal D* **68**, 272 (2014), doi: `10.1140/epjd/e2014-50500-1 <http://dx.doi.org/10.1140/epjd/e2014-50500-1>`_

