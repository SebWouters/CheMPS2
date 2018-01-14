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

.. index:: Symmetry

Symmetries
==========

CheMPS2 exploits the :math:`\mathsf{SU(2)}` spin symmetry, :math:`\mathsf{U(1)}` particle number symmetry, and abelian point group symmetries :math:`\{ \mathsf{C_1}, \mathsf{C_i}, \mathsf{C_2}, \mathsf{C_s}, \mathsf{D_2}, \mathsf{C_{2v}}, \mathsf{C_{2h}}, \mathsf{D_{2h}}  \}` of ab initio quantum chemistry Hamiltonians. Thereto the orbital occupation and virtual indices have to be represented by states which transform according to a particular row of one of the irreps of the symmetry group of the Hamiltonian. For example for orbital :math:`k`:

.. math::

    \left|-\right\rangle & \rightarrow & \left|s = 0;s^z=0;N=0; I=I_0\right\rangle\\
    \left|\uparrow\right\rangle & \rightarrow & \left|s = \frac{1}{2};s^z=\frac{1}{2};N=1; I=I_k\right\rangle\\
    \left|\downarrow\right\rangle & \rightarrow & \left|s = \frac{1}{2};s^z=-\frac{1}{2};N=1; I=I_k\right\rangle\\
    \left|\uparrow\downarrow\right\rangle & \rightarrow & \left|s = 0;s^z=0;N=2; I=I_k \otimes I_k = I_0\right\rangle.

Then the MPS tensors factorize into Clebsch-Gordan coefficients and reduced tensors due to the Wigner-Eckart theorem:

.. math::

    A[i]^{(ss^zNI)}_{(j_L j_L^z N_L I_L \alpha_L);(j_R j_R^z N_R I_R \alpha_R)} = \left\langle j_L j_L^z s s^z \mid j_R j_R^z \right\rangle \delta_{N_L+N,N_R} \delta_{I_L\otimes I, I_R} T[i]^{(sNI)}_{(j_L N_L I_L \alpha_L);(j_R N_R I_R \alpha_R)}.

This has three important consequences:

#. There is block-sparsity due to the Clebsch-Gordan coefficients. Remember that the Clebsch-Gordan coefficients of abelian groups are Kronecker :math:`\delta`'s. The block-sparsity results in both memory and CPU time savings.
#. There is information compression for spin symmetry sectors other than singlets, as the tensor :math:`\mathbf{T[i]}` does not contain spin projection indices. The virtual dimension associated with :math:`\mathbf{T[i]}` is called the reduced virtual dimension :math:`D_{\mathsf{SU(2)}}`. This also results in both memory and CPU time savings.
#. Excited states in different symmetry sectors can be obtained by ground-state calculations.

The operators

.. math::

    \hat{b}^{\dagger}_{k\sigma} & = & \hat{a}^{\dagger}_{k\sigma}\\
    \hat{b}_{k\sigma} & = & (-1)^{\frac{1}{2}-\sigma}\hat{a}_{k-\sigma}

of orbital :math:`c` transform according to row :math:`(s = \frac{1}{2}; s^z=\sigma; N=\pm 1; I_c)` of irrep :math:`(s = \frac{1}{2}; N=\pm 1; I_c)`. :math:`\hat{b}^{\dagger}` and :math:`\hat{b}` are hence both doublet irreducible tensor operators, and the Wigner-Eckart theorem allows to factorize corresponding matrix elements into Clebsch-Gordan coefficients and reduced matrix elements. Together with the Wigner-Eckart theorem for the MPS tensors, this allows to work with reduced quantities only in CheMPS2. Only Wigner 6-j and 9-j symbols are needed, but never Wigner 3-j symbols or Clebsch-Gordan coefficients.

For more information on the exploitation of symmetry in the DMRG method, please read Ref. [SYMM1]_.

.. [SYMM1] S. Wouters and D. Van Neck, *European Physical Journal D* **68**, 272 (2014), doi: `10.1140/epjd/e2014-50500-1 <http://dx.doi.org/10.1140/epjd/e2014-50500-1>`_

