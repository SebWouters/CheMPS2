.. index:: Computer resources
.. index:: Wall time
.. index:: Memory
.. index:: Disk

Typical resource requirements
=============================

In this section, typical resource requirements for DMRG calculations are discussed. With :math:`L` spatial orbitals and :math:`D` virtual basis states, the algorithm has a theoretical scaling per sweep of

* :math:`\mathcal{O}(L^4D^2 + L^3D^3)` in CPU time
* :math:`\mathcal{O}(L^2D^2)` in memory
* :math:`\mathcal{O}(L^3D^2)` in disk

The block-sparsity and information compression due to the exploitation of symmetry have not been taken into account in these scalings! Ref. [TIMING]_ contains CPU time measurements for polyenes of increasing length, and demonstrates the scaling of CheMPS2 with :math:`L`.


N2/cc-pVDZ
----------

The nitrogen dimer in the cc-pVDZ basis has an active space of 14 electrons in 28 orbitals. The exploited point group in the calculations was :math:`\mathsf{d2h}`, and the targeted state was :math:`X^1\Sigma_g^+` at equilibrium bond length: 2.118 a.u. This system was first studied with DMRG in Ref. [NITROGEN]_. The listed CheMPS2 timings are wall times per sweep (in seconds) on 16 Intel Xeon Sandy Bridge (E5-2670) cores @ 2.6 GHz. The calculation was performed on a single node, and needed ~ 6 Gb of memory.

 +----------------------------+-------------------------+--------------------+-----------------------+
 | :math:`D_{\mathsf{SU(2)}}` | Wall time per sweep (s) | :math:`w_D^{disc}` | :math:`E_D` (Hartree) |
 +============================+=========================+====================+=======================+
 | 1000                       | 117                     | 9.9532e-07         | -109.28209681         |
 +----------------------------+-------------------------+--------------------+-----------------------+
 | 1500                       | 279                     | 3.7204e-07         | -109.28214478         |
 +----------------------------+-------------------------+--------------------+-----------------------+
 | 2000                       | 529                     | 1.7795e-07         | -109.28216006         |
 +----------------------------+-------------------------+--------------------+-----------------------+
 | 2500                       | 1339                    | 9.6998e-08         | -109.28216623         |
 +----------------------------+-------------------------+--------------------+-----------------------+


H2O/Roos' ANO DZ
----------------

Water in Roos' ANO DZ basis has an active space of 10 electrons in 41 orbitals. The exploited point group in the calculations was :math:`\mathsf{c2v}`, and the targeted state was :math:`^1A_1` at equilibrium geometry: O @ (0, 0, 0) and H @ (± 0.790689766, 0, 0.612217330) Angstrom. This system was first studied with DMRG in Ref. [WATER]_. The listed CheMPS2 timings are wall times per sweep (in seconds) on 20 Intel Xeon Ivy Bridge (E5-2670 v2) cores @ 2.5 GHz. The calculation was performed on a single node, and needed ~ 64 Gb of memory.

 +----------------------------+-------------------------+--------------------+-----------------------+
 | :math:`D_{\mathsf{SU(2)}}` | Wall time per sweep (s) | :math:`w_D^{disc}` | :math:`E_D` (Hartree) |
 +============================+=========================+====================+=======================+
 | 1000                       | 820                     | 8.6093e-08         | -76.31468322          |
 +----------------------------+-------------------------+--------------------+-----------------------+
 | 2000                       | 4304                    | 1.1262e-08         | -76.31471043          |
 +----------------------------+-------------------------+--------------------+-----------------------+
 | 3000                       | 11872                   | 2.8027e-09         | -76.31471338          |
 +----------------------------+-------------------------+--------------------+-----------------------+
 | 4000                       | 23915                   | 6.9927e-10         | -76.31471401          |
 +----------------------------+-------------------------+--------------------+-----------------------+


.. [TIMING] S. Wouters and D. Van Neck, *European Physical Journal D* **68**, 272 (2014), doi: `10.1140/epjd/e2014-50500-1 <http://dx.doi.org/10.1140/epjd/e2014-50500-1>`_
.. [NITROGEN] G.K.-L. Chan, M. Kallay and J. Gauss, *Journal of Chemical Physics* **121**, 6110 (2004), doi: `10.1063/1.1783212 <http://dx.doi.org/10.1063/1.1783212>`_
.. [WATER] G. K.-L. Chan and M. Head-Gordon, *Journal of Chemical Physics* **118**, 8551 (2003), doi: `10.1063/1.1574318 <http://dx.doi.org/10.1063/1.1574318>`_


