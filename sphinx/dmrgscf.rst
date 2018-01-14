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

.. index:: DMRG-SCF
.. index:: DIIS
.. index:: Newton-Raphson
.. index:: augmented Hessian

DMRG-SCF
========

In methods which use a FCI solver, this solver can be replaced by DMRG. DMRG allows for an efficient extraction of the 2-RDM. The 2-RDM of the active space is required in the complete active space self-consistent field (CASSCF) method to compute the gradient and the Hessian with respect to orbital rotations [DMRGSCF1]_. It is therefore natural to introduce a CASSCF variant with DMRG as active space solver, called DMRG-SCF [DMRGSCF2]_ [DMRGSCF3]_ [DMRGSCF4]_, which allows to treat static correlation in large active spaces. In CheMPS2, the augmented Hessian Newton-Raphson DMRG-SCF method is implemented, with exact Hessian [DMRGSCF5]_ [DMRGSCF6]_.

Augmented Hessian Newton-Raphson
--------------------------------
    
The basic idea is to express the energy with the unitary group generators up to second order:

.. math::

    \hat{E}_{pq} & = & \sum\limits_{\sigma} \hat{a}^{\dagger}_{p \sigma} \hat{a}_{q \sigma} \\
    \left[ \hat{E}_{pq} , \hat{E}_{rs} \right] & = & \delta_{qr} \hat{E}_{ps} - \delta_{ps} \hat{E}_{rq} \\
    \hat{E}^{-}_{pq} & = & \hat{E}_{pq} - \hat{E}_{qp} \\
    \hat{T}(\vec{x}) & = & \sum\limits_{p<q} x_{pq} \hat{E}^{-}_{pq} \\
    E(\vec{x}) & = & \left\langle 0 \mid e^{\hat{T}(\vec{x})} \hat{H} e^{-\hat{T}(\vec{x})} \mid 0 \right\rangle \approx E(0) + \vec{x}^T \vec{g} + \frac{1}{2} \vec{x}^T \mathbf{H} \vec{x}

The vector :math:`\vec{g}` is the gradient and the matrix :math:`\mathbf{H}` the Hessian for orbital rotations. The minimum of :math:`E(\vec{x})` is found at :math:`\vec{x} = - \mathbf{H}^{-1} \vec{g}`. The variables :math:`\vec{x}` parametrize an additional orbital rotation :math:`\mathbf{U}_{add} = \exp(\mathbf{T}(\vec{x}))`, with :math:`\mathbf{T}(\vec{x}) = -\mathbf{T}^T(\vec{x})` a real-valued skew-symmetric matrix. The additional orbital rotation :math:`\mathbf{U}_{add}` transforms the current orbitals :math:`\mathbf{U}(n)` to the new orbitals

.. math::

    \mathbf{U}(n+1) = \mathbf{U}_{add} \mathbf{U}(n) = \exp(\mathbf{T}(\vec{x}(n))) \mathbf{U}(n).

This updating scheme is called the Newton-Raphson method. If the Hessian is positive definite, these updates are stable. For saddle points in the energy landscape, the Hessian has negative eigenvalues, and these updates can be unstable. It is therefore better to use the augmented Hessian Newton-Raphson method:

.. math::

    \left[ \begin{array}{cc} \mathbf{H} & \vec{g} \\ \vec{g}^T & 0 \end{array} \right] \left[ \begin{array}{c} \vec{x} \\ 1 \end{array} \right] = \alpha \left[ \begin{array}{c} \vec{x} \\ 1 \end{array} \right].

The eigenvector with smallest algebraic eigenvalue determines a stable update :math:`\vec{x}`, as is well explained in Ref. [DMRGSCF7]_.

DIIS
----

When the update norm :math:`\|\vec{x}\|_2` is small enough, the convergence can be accelerated by the direct inversion of the iterative subspace (DIIS) [DMRGSCF8]_ [DMRGSCF9]_. For a given set of orbitals :math:`\mathbf{U}(n)`, the update :math:`\vec{x}(n)` is calculated with the augmented Hessian Newton-Raphson method. This update defines the next set of orbitals:

.. math::

    \mathbf{U}(n+1) = \mathbf{U}_{add} \mathbf{U}(n) = \exp(\mathbf{T}(\vec{x}(n))) \mathbf{U}(n).

In DIIS, the error vector :math:`\vec{x}(n)` and the state vector :math:`\mathbf{Y}(n+1) = \log(\mathbf{U}(n+1))` are added to a list. The error

.. math::

    e = \left\| \sum\limits_{i = 1}^{\kappa} c_i \vec{x}(n - \kappa + i) \right\|_2

is minimized under the constraint :math:`\sum_{i} c_i = 1`. :math:`\kappa` is the size of the list memory, i.e. the number of retained vectors. The minimization of the error :math:`e` can be performed with Lagrangian calculus:

.. math::

    \left[ \begin{array}{cc} \mathbf{B} & \vec{1} \\ \vec{1}^T & 0 \end{array} \right] \left[ \begin{array}{c} \vec{c} \\ \lambda \end{array} \right] = \left[ \begin{array}{c} \vec{0} \\ 1 \end{array} \right],

where :math:`2\lambda` is the Lagrangian multiplier and

.. math::

    \left[\mathbf{B}\right]_{ij} = \vec{x}^T(n - \kappa + i) \vec{x}(n - \kappa + j) = \left[\mathbf{B}\right]_{ji}.

The new state vector is then defined as

.. math::

    \mathbf{Y}_{new} = \sum\limits_{i = 1}^{\kappa} c_i \mathbf{Y}(n+1 - \kappa + i).

The current orbitals are then set to :math:`\mathbf{U}(n+1) = \exp(\mathbf{Y}_{new})`.

.. [DMRGSCF1] P.E.M. Siegbahn, J. Almlof, A. Heiberg and B.O. Roos, *Journal of Chemical Physics* **74**, 2384-2396 (1981), doi: `10.1063/1.441359 <http://dx.doi.org/10.1063/1.441359>`_
.. [DMRGSCF2] D. Ghosh, J. Hachmann, T. Yanai and G.K.-L. Chan, *Journal of Chemical Physics* **128**, 144117 (2008), doi: `10.1063/1.2883976 <http://dx.doi.org/10.1063/1.2883976>`_
.. [DMRGSCF3] D. Zgid and M. Nooijen, *Journal of Chemical Physics* **128**, 144116 (2008), doi: `10.1063/1.2883981 <http://dx.doi.org/10.1063/1.2883981>`_
.. [DMRGSCF4] T. Yanai, Y. Kurashige, D. Ghosh and G.K.-L. Chan, *International Journal of Quantum Chemistry* **109**, 2178-2190 (2009), doi: `10.1002/qua.22099 <http://dx.doi.org/10.1002/qua.22099>`_
.. [DMRGSCF5] S. Wouters, W. Poelmans, P.W. Ayers and D. Van Neck, *Computer Physics Communications* **185**, 1501-1514 (2014), doi: `10.1016/j.cpc.2014.01.019 <http://dx.doi.org/10.1016/j.cpc.2014.01.019>`_
.. [DMRGSCF6] S. Wouters, T. Bogaerts, P. Van Der Voort, V. Van Speybroeck and D. Van Neck, *Journal of Chemical Physics* **140**, 241103 (2014), doi: `10.1063/1.4885815 <http://dx.doi.org/10.1063/1.4885815>`_
.. [DMRGSCF7] A. Banerjee, N. Adams, J. Simons and R. Shepard, *Journal of Physical Chemistry* **89**, 52-57 (1985), doi: `10.1021/j100247a015 <http://dx.doi.org/10.1021/j100247a015>`_
.. [DMRGSCF8] P. Pulay, *Chemical Physics Letters* **73**, 393-398 (1980), doi: `10.1016/0009-2614(80)80396-4 <http://dx.doi.org/10.1016/0009-2614(80)80396-4>`_
.. [DMRGSCF9] C.D. Sherrill, Programming DIIS, http://vergil.chemistry.gatech.edu/notes/diis/node3.html (2000).

