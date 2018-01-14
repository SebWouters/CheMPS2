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

.. index:: CASPT2

Internally contracted CASPT2
============================

The required 4-RDM elements
---------------------------

As DMRG is a wavefunction method, higher order reduced density matrices can be evaluated.
In CheMPS2, up to the 3-RDM is implemented with the usual sweep algorithm.
For internally contracted CASPT2, the contraction of the 4-RDM with the generalized Fock operator is needed. The generalized Fock operator

.. math::

    \hat{F} & = & \sum\limits_{pq} F_{pq} \hat{E}_{pq}

has matrix elements

.. math::

    F_{pq}  =  \frac{1}{2} \sum\limits_{\sigma} \left\langle \hat{a}_{p\sigma} \left[ \hat{H}, \hat{a}_{q \sigma}^{\dagger} \right] - \hat{a}_{p\sigma}^{\dagger} \left[ \hat{H}, \hat{a}_{q \sigma} \right]  \right\rangle  =  t_{pq} + \sum\limits_{rs} \left\langle \hat{E}_{rs} \right\rangle \left( \left( pq | rs \right) - \frac{1}{2} \left( pr | qs \right) \right).

These matrix elements are symmetric and diagonal in the spatial orbital irreps. The required contraction of the 4-RDM with the generalized Fock operator

.. math::

     \left( F . \Gamma^4 \right)_{pqr; wxy} & = & \sum_{\sigma \tau \upsilon} \left\langle \hat{a}^{\dagger}_{p \sigma} \hat{a}^{\dagger}_{q \tau} \hat{a}^{\dagger}_{r \upsilon} \hat{F} \hat{a}_{y \upsilon} \hat{a}_{ x \tau} \hat{a}_{w \sigma} \right\rangle

can therefore be obtained from excited wavefunctions, formed by symmetry-conserving single-particle excitations op top of the reference wavefunction. For the matrix element :math:`F_{sz} = F_{zs}` the following sum of 4-RDM elements is needed:

.. math::

    \Gamma^4_{pqrs; wxyz} + \Gamma^4_{pqrz; wxys}.

It can be obtained by calculating

.. math::

    \left| sz, \alpha, \beta \right\rangle = \left[ \alpha \left( \hat{E}_{sz} + {E}_{zs} \right) + \beta \right]  \left| \Psi_0 \right\rangle

in which the spatial orbitals :math:`s` and :math:`z` have equal irreps. This excited wavefunction is decomposed into an MPS, with a sweep algorithm with negligible cost [CASPT2]_. Denote the 3-RDM of the (unnormalized) excited wavefunctions as

.. math::

    \Gamma(sz,\alpha,\beta)^3_{pqr; wxy} = \sum\limits_{ \sigma \tau \upsilon } \left\langle sz, \alpha, \beta \mid \hat{a}_{p\sigma}^{\dagger} \hat{a}_{q\tau}^{\dagger} \hat{a}_{r\upsilon}^{\dagger} \hat{a}_{y\upsilon} \hat{a}_{x\tau} \hat{a}_{w\sigma} \mid sz, \alpha, \beta \right\rangle.

In this notation :math:`\Gamma( sz, 0, 1 )^3_{pqr; wxy} = \Gamma^3_{pqr; wxy}`, the 3-RDM of the reference wavefunction. The following identity holds:

.. math::

    & & 2 \left[ \Gamma^4_{pqrs; wxyz} + \Gamma^4_{pqrz; wxys} \right] + \Gamma^3_{pqr; wxy} \\
    & = & \Gamma( sz, 1, 1 )^3_{pqr; wxy} - \Gamma( sz, 1, 0 )^3_{pqr; wxy} \\
    & - & \delta_{s,p} \Gamma^3_{zqr; wxy} - \delta_{s,q} \Gamma^3_{pzr; wxy} - \delta_{s,r} \Gamma^3_{pqz; wxy} \\
    & - & \delta_{s,w} \Gamma^3_{pqr; zxy} - \delta_{s,x} \Gamma^3_{pqr; wzy} - \delta_{s,y} \Gamma^3_{pqr; wxz} \\
    & - & \delta_{z,p} \Gamma^3_{sqr; wxy} - \delta_{z,q} \Gamma^3_{psr; wxy} - \delta_{z,r} \Gamma^3_{pqs; wxy} \\
    & - & \delta_{z,w} \Gamma^3_{pqr; sxy} - \delta_{z,x} \Gamma^3_{pqr; wsy} - \delta_{z,y} \Gamma^3_{pqr; wxs},

which allows to obtain the contraction :math:`\left( F . \Gamma^4 \right)` by calculating the 3-RDM of several excited wavefunctions! This algorithm has in practice the same computational cost as the regular 4-RDM evaluation during the usual sweep algorithm.

CAS perturbation theory
-----------------------

The full Hilbert space :math:`\mathcal{H}` is split up into four parts [ROOS1]_ [ROOS2]_:

.. math::
    \mathcal{H} = \mathcal{V}_0 \oplus \mathcal{V}_{\text{K}} \oplus \mathcal{V}_{\text{SD}} \oplus \mathcal{V}_{\text{TQ..}}.

#. :math:`\mathcal{V}_0` contains only the CASSCF solution :math:`\left| \Psi_0 \right\rangle`.
#. :math:`\mathcal{V}_{\text{K}}` is the space spanned by all possible active space excitations on top of :math:`\left| \Psi_0 \right\rangle` which are orthogonal to :math:`\mathcal{V}_0`. Wavefunctions in :math:`\mathcal{V}_{\text{K}}` have the same core and virtual orbitals as :math:`\left| \Psi_0 \right\rangle`, with the same occupation.
#. :math:`\mathcal{V}_{{\text{SD}}}` contains all single and double particle excitations on top of :math:`\left| \Psi_0 \right\rangle` which are orthogonal to :math:`\mathcal{V}_0 \oplus \mathcal{V}_{\text{K}}`. With the indices :math:`ij` for core orbitals, :math:`tuv` for active orbitals, and :math:`ab` for virtual orbitals, :math:`\mathcal{V}_{{\text{SD}}}` is spanned by the following excitation types:

    .. math::

        \text{A} & : & \quad \hat{E}_{ti} \hat{E}_{uv} \left| \Psi_0 \right\rangle, \\
        \text{B} & : & \quad \hat{E}_{ti} \hat{E}_{uj} \left| \Psi_0 \right\rangle, \\
        \text{C} & : & \quad \hat{E}_{at} \hat{E}_{uv} \left| \Psi_0 \right\rangle, \\
        \text{D} & : & \quad \hat{E}_{ai} \hat{E}_{tu} \left| \Psi_0 \right\rangle,~\hat{E}_{ti}\hat{E}_{au} \left| \Psi_0 \right\rangle, \\
        \text{E} & : & \quad \hat{E}_{ti} \hat{E}_{aj} \left| \Psi_0 \right\rangle, \\
        \text{F} & : & \quad \hat{E}_{at} \hat{E}_{bu} \left| \Psi_0 \right\rangle, \\
        \text{G} & : & \quad \hat{E}_{ai} \hat{E}_{bt} \left| \Psi_0 \right\rangle, \\
        \text{H} & : & \quad \hat{E}_{ai} \hat{E}_{bj} \left| \Psi_0 \right\rangle.

#. And :math:`\mathcal{V}_{{\text{TQ..}}}` is the remainder of :math:`\mathcal{H}`.

The zeroth order Hamiltonian for internally contracted CASPT2 is

.. math::
    \hat{H}_0 = \hat{P}_0 \hat{F} \hat{P}_0 + \hat{P}_{\text{K}} \hat{F} \hat{P}_{\text{K}} + \hat{P}_{{\text{SD}}} \hat{F} \hat{P}_{{\text{SD}}} + \hat{P}_{{\text{TQ..}}} \hat{F} \hat{P}_{{\text{TQ..}}},

where :math:`\hat{P}_{\text{X}}` is the projector onto :math:`\mathcal{V}_{\text{X}}`.
The first order wavefunction :math:`\left| \Psi_1 \right\rangle` for internally contracted CASPT2 is spanned by a linear combination over :math:`\mathcal{V}_{\text{SD}}`:

.. math::
    \left| \Psi_1 \right\rangle = \sum_{pq;rs \in \mathcal{V}_{\text{SD}}} C_{pq;rs} \hat{E}_{pq} \hat{E}_{rs} \left| \Psi_0 \right\rangle = \sum_{pq;rs \in \mathcal{V}_{\text{SD}}} C_{pq;rs} \left| \Psi_{pq;rs} \right\rangle.

The coefficients can be found by solving

.. math::
    \sum_{pq;rs \in \mathcal{V}_{\text{SD}}} \left\langle \Psi_{wx;yz} \mid \hat{H}_0 - E_0 \mid \Psi_{pq;rs} \right\rangle C_{pq;rs} = - \left\langle \Psi_{wx;yz} \mid \hat{H} \mid \Psi_0 \right\rangle.

The overlap matrix :math:`\left\langle \Psi_{wx;yz} \mid \Psi_{pq;rs} \right\rangle` is block-diagonal in the different excitation types (A to H). It is diagonalized, small eigenvalues are discarded, and the linear equation is transformed to

.. math::

    \sum\limits_{ \beta } \left( \mathcal{F}_{\alpha\beta} - E_0 \delta_{\alpha,\beta} \right) \mathcal{C}_{\beta} = - \mathcal{V}_{\alpha},

with :math:`\mathcal{F}` diagonal for two excitations of the same type (A to H). The following initial guess is used to solve this linear equation with either the conjugate gradient or Davidson algorithm:

.. math::

    \mathcal{C}_{\alpha}^{\text{ini}} = - \frac{ \mathcal{V}_{\alpha} }{ \mathcal{F}_{\alpha\alpha} - E_0 }.

If the active space orbitals in the DMRG algorithm are not pseudocanonical, :math:`\Gamma^1`, :math:`\Gamma^2`, :math:`\Gamma^3`, and :math:`\left(F.\Gamma^4\right)` are rotated to the pseudocanonical orbital basis before building the required intermediates to solve the CASPT2 linear equation.

In order to mitigate intruder state problems, CheMPS2 allows to specify an imaginary level shift [IMAG]_ and/or ionization potential - electron affinity shift [IPEA]_. For the latter, the left-hand side matrix of the CASPT2 linear equation is shifted with

.. math::

   \left\langle \Psi_{ wx;yz } \mid \hat{F} \mid \Psi_{pq;rs } \right\rangle \mathrel{+}= \delta_{p,w} \delta_{q,x} \delta_{r,y} \delta_{s,z} \frac{ \epsilon^{\text{IPEA}}}{2} \left\langle \Psi_{wx;yz} \mid \Psi_{pq;rs} \right\rangle \left( 4 + \left\langle \hat{E}_{pp} \right\rangle - \left\langle \hat{E}_{qq} \right\rangle + \left\langle \hat{E}_{rr} \right\rangle - \left\langle \hat{E}_{ss} \right\rangle \right).


CASPT2 calculations
-------------------

In order to calculate the CASPT2 variational second order perturbation correction energy, the following call should be made:

.. code-block:: c++

    double CheMPS2::CASSCF::caspt2( const int Nelectrons, const int TwoS, const int Irrep, ConvergenceScheme * OptScheme, const int rootNum, DMRGSCFoptions * scf_options, const double IPEA, const double IMAG, const bool PSEUDOCANONICAL, const bool CHECKPOINT, const bool CUMULANT )

**after** the CASSCF orbital rotations are converged.

#. The first six parameters are the same as for ``CheMPS2::CASSCF::solve`` in :ref:`CheMPS2::CASSCF <label-casscf-calculations-api>`.
#. IPEA and IMAG are the ionization potential - electron affinity and imaginary level shifts.
#. PSEUDOCANONICAL allows to change the converged CASSCF orbitals (localized, natural, ...) to pseudocanonical orbitals **before** the DMRG calculation to obtain the contraction of the 4-RDM with the generalized Fock operator. This has the advantage that the Fock operator matrix elements are diagonal, which leads to a significant reduction in computational cost **if** it not requires a significantly larger virtual dimension :math:`D`.
#. CHECKPOINT allows to switch on the creation of checkpoints to calculate the required contraction during multiple runs. If CHECKPOINT is true, then after the initial run

    .. code-block:: c++

        scf_options->setDoDIIS( false )
        scf_options->setWhichActiveSpace( 0 )

   should be set in order to use exactly the same orbitals in consecutive runs!
#. It is advised to leave CUMULANT = false.


.. [CASPT2] S. Wouters, V. Van Speybroeck and D. Van Neck, *Journal of Chemical Physics* **145**, 054120 (2016), doi: `10.1063/1.4959817 <http://dx.doi.org/10.1063/1.4959817>`_
.. [ROOS1]  K. Andersson, P.-A. Malmqvist, B.O. Roos, A.J. Sadlej and K. Wolinski, *Journal of Physical Chemistry* **94**, 5483-5488 (1990). doi: `10.1021/j100377a012 <http://dx.doi.org/10.1021/j100377a012>`_
.. [ROOS2]  K. Andersson, P.‐A. Malmqvist and B.O. Roos, *Journal of Chemical Physics* **96**, 1218-1226 (1992). doi: `10.1063/1.462209 <http://dx.doi.org/10.1063/1.462209>`_
.. [IMAG]   N. Forsberg and P.-A. Malmqvist, *Chemical Physics Letters* **274**, 196-204 (1997). doi: `10.1016/S0009-2614(97)00669-6 <http://dx.doi.org/10.1016/S0009-2614(97)00669-6>`_
.. [IPEA]   G. Ghigo, B.O. Roos and P.-A. Malmqvist, *Chemical Physics Letters* **396**, 142-149 (2004). doi: `10.1016/j.cplett.2004.08.032 <http://dx.doi.org/10.1016/j.cplett.2004.08.032>`_


