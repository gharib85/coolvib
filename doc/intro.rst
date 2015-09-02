Introduction
============

What is the purpose of *coolvib*?
---------------------------------

**coolvib** is a `Python <www.python.org>`_-based package that enables the 
calculation of vibrational cooling and friction due to non-adiabatic coupling 
between electrons and vibrations or phonons. The electronic structure is thereby 
treated as a perturbation acting on the nuclear motion. **coolvib** is written 
as a post-processing tool that interfaces with different codes (see :ref:`parser`).

**coolvib** facilitates calculation of input data using these codes and heavily depends 
on interfaces and data structures as provided by the `Atomic Simulation Environment <https://wiki.fysik.dtu.dk/ase/>`_.

What is vibrational cooling and when is it important?
--------------------------------------------------------

*Vibrational cooling* [Wodtke2004]_, [Gergen2001]_ in general refers to the fact that vibrations or nuclear motion of atoms and molecules can be affected by external perturbations leading to energy exchange and damping or cooling down of these affected motions. Such effects are specifically relevant for atoms and molecules in contact with many degrees of freedom such as when in contact with surfaces or nanostructures. The two dominantly discussed reasons for vibrational cooling or energy dissipation are related to coupling of vibrations and phonons with each other and coupling of vibrations with electronic degrees of freedom. **coolvib** targets the latter, namely the **non-adiabatic vibrational cooling and frictional force** acting on adsorbate atoms and molecules due to coupling between electrons and nuclei [Head-Gordon1992]_, [Head-Gordon1995]_, [Persson1982]_. 

Non-adiabatic friction is specifically important for atoms and molecules adsorbed at metal surfaces. In such systems nuclear motion and electronic excitations can be induced by energies of similar scale leading to a non-vanishing coupling between the two. This is equivalent to a break-down of the Born-Oppenheimer approximation. [Wodtke2004]_


How does *coolvib* treat non-adiabatic friction and vibrational cooling?
------------------------------------------------------------------------

**coolvib** uses Time Dependent Perturbation Theory [Persson1982]_ [Head-Gordon1992]_ to calculate the tensorial friction acting on adsorbate motion. It thereby assumes that the effect that electronic excitations have on the nuclear motion is *small* and can therefore be treated as a perturbation. The result of this calculation is a decay rate along a given vibrational mode or a full cartesian friction tensor from which vibrational lifetimes can be derived.

For more detail on how this is calculated see the  :ref:`theory` section of this documentation.


References
----------

.. [Wodtke2004] `A.` Wodtke, J. Tully, and D. Auerbach, *Int. Rev. Phys. Chem.* **23**, 513-539 (2004)
.. [Gergen2001] `B.` Gergen, H. Nienhaus, W. H. Weinberg, and E. W. McFarland, *Science* **294**, 2521-2523 (2001)
.. [Head-Gordon1992] `M.` Head-Gordon, and J. Tully, *J. Chem. Phys* **96**, 3939 (1992)
.. [Head-Gordon1995] `M.` Head-Gordon, and J. Tully, *J. Chem. Phys* **103**, 10137 (1995)
.. [Persson1982] `M.` Persson, and B. Hellsing, *Phys. Rev. Lett.* **49**, 662-665 (1982)
