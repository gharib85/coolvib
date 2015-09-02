.. _derivaton:

Derivation
==========

The following derivations are directly taken from [Maurer2015]_.

The electron-phonon interaction
-------------------------------

Our starting point is the molecular many-body problem in the 
Born-Oppenheimer approximation. 
Herein we define an electronic Schrödinger equation (SGE):

.. math::
     \hat{H}_{\mathrm{e}}\Phi(t,\mathbf{r};\mathbf{R}) = i\hbar\frac{d}{dt} \Phi(t,\mathbf{r};\mathbf{R})

with

.. math::
   :label: Helec
  
    \hat{H}_{\mathrm{e}} &= -\sum_{i=1}^N\frac{\hbar^2}{2}\nabla_i^2 +\sum_{i=1}^N\sum_{j>i}^N\frac{1}{r_{ij}}-\sum_{i=1}^N\sum_{A=1}^M \frac{Z_A}{r_{iA}}+ \\ \nonumber &+ \sum_{A=1}^M\sum_{B>A}^M \frac{Z_AZ_B}{R_{AB}}

and a nuclear Schrödinger equation as

.. math::
   \hat{H}_{\mathrm{n}}\chi(t,\mathbf{R}) = i\hbar\frac{d}{dt} \chi(t,\mathbf{R})

with

.. math::
    
   \hat{H}_{\mathrm{n}} = -\sum_{A=1}^M\frac{\hbar^2}{2M_A}\nabla_A^2  + E_{\mathrm{e}}(\mathbf{R}) .
    
All of the above is given in atomic units, with the exception of explicitly stating :math:`\hbar(=1)`.
The electronic SGE and its eigenfunctions depend parametrically on the nuclear coordinates. 
Due to separation of time-scales at which electrons and nuclei move, these equations are solved independently.

We can introduce coupling between the two systems by constructing the full molecular wavefunction 
as linear combination of products of nuclear and adiabatic electronic wavefunctions in a 
so-called Born expansion

.. math::
    \Psi(t,\mathbf{r},\mathbf{R}) = \sum_j \chi_j(t,\mathbf{R}) \Phi_j(\mathbf{r};\mathbf{R}) 

The coupling is introduced by the third term in eq. :eq:`Helec`. 
Insertion into the full time-dependent Schrödinger equation, 
multiplication from the left with an electronic wavefunction :math:`\Phi_i` 
and integration over electronic degrees of freedom we find an additional term in the nuclear Schrödinger equation:
    
.. math::
     (\underbrace{\hat{T}_{{\mathrm{n}}} + E^{\mathrm{e}}_j}_{\hat{H}_{\mathrm{n}}}+ \hat{H}_{\mathrm{ep}})\chi_j = i\hbar\frac{d}{dt}\chi_j

with the action of the electron-phonon coupling term on the nuclear wavefunction defined as

.. math::
   :label: app-Hep
    
    \hat{H}_{\mathrm{ep}}\chi_j =-\frac{\hbar^2}{2M}\sum_{i\neq j}\left[2\langle\Phi_i|\nabla_{\mathbf{R}}\Phi_j\rangle\cdot\nabla_{\mathbf{R}}+\langle\Phi_i|\nabla^2_{\mathbf{R}}\Phi_j\rangle \right]\chi_i
    

One way to understand this interaction is by viewing it as hybridization between 
electronic and nuclear wavefunctions. 
Each nuclear wavefunction is split into as many states as there are electronic excitations (and *vice versa*).
Throughout the remainder of this documentation we will assume :math:`\hat{H}_{\mathrm{ep}}` to be small 
compared to :math:`\hat{H}_n`, i.e. it does not change the nuclear dynamics qualitatively. 
This enables application of perturbation theory. In the following we will find an expression for the 
interaction matrix elements and subsequently use Time Dependent Perturbation Theory
to derive our working expressions.


Derivation of coupling matrix elements
--------------------------------------

Next we will find an expression for the matrix elements

.. math::

    g_{fj} = \braket{\chi_f|\hat{H}_{\mathrm{ep}}|\chi_j} = -\frac{\hbar^2}{2M}\sum_{i}\bra{\chi_f}\biggl[2\braket{\Phi_i|\nabla_{\mathbf{R}}\Phi_j}\cdot\nabla_{\mathbf{R}}  +\braket{\Phi_i|\nabla^2_{\mathbf{R}}\Phi_j} \biggr]\ket{\chi_i}.

Our starting point is a periodic condensed phase material in which the electronic structure underlying the nuclear motion is given by effectively independent quasi-particles such as calculated with semi-local approximations to Density Functional Theory or the GW method. The unperturbed nuclear motion is described in terms of harmonic vibrations as given by the following Hamiltonian:

.. math::
   :label: app-HO
    
    H_0 = \sum_{\mathbf{q}j} \hbar \omega_{\mathbf{q}j}\left(b^\dagger_{\mathbf{q}j}b_{\mathbf{qj}} +\frac{1}{2} \right) .

Herein :math:`b^\dagger_{\mathbf{q}j}` and :math:`b_{\mathbf{qj}}` are the usual particle annihilation and creation operators for a phonon :math:`j` at wavevector :math:`\mathbf{q}`. Phonon normal mode displacement vectors are given as :math:`\mathbf{e}_{\mathbf{q}j}`. With this choice the scalar coupling term in eq. :eq:`app-Hep` vanishes due to orthonormality of the nuclear wavefunctions. The matrix element between two many-body nuclear wavefunction states, which only differ by a single one-particle state reduces to 

.. math::
    
    \braket{\chi_f|\hat{H}_{\mathrm{ep}}|\chi_j} \Rightarrow \braket{\chi_{\mathbf{q}f}|\hat{H}_{\mathrm{ep}}|\chi_{\mathbf{q}j}}.

Furthermore, by taking the many-body electronic wavefunction as a product of quasi-particle states with quantum states specified by the quantum numbers :math:`s`, :math:`\mathbf{k}`, and :math:`\nu` representing spin, wave vector, and band index, all possible single-particle excited states are defined with the above indices for electrons and phonons. Correspondingly the matrix elements are

.. math:: 
    \braket{\chi_{\mathbf{q}f}|\hat{H}_{\mathrm{ep}}|\chi_{\mathbf{q}j}} = -\frac{\hbar^2}{M} \sum_{i}\cdot \braket{\chi_{\mathbf{q}f}|\sum_{s}\sum_{\mathbf{k}\nu\mathbf{k}'\nu'}  \braket{\psi_{s\mathbf{k}\nu}|\mathbf{e}_{\mathbf{q}j}\nabla_{\mathbf{R}}\psi_{s\mathbf{k}'\nu'}}\mathbf{e}_{\mathbf{q}j}\nabla_{\mathbf{R}}|\chi_{\mathbf{q}i}} ,

where the different sums run over occupied and unoccupied states in the ground state, respectively. The occupations of the one-particle electronic states and thereby the Fermi level are defined by Fermi-Dirac statistics. With harmonic displacements as starting point we have defined the nuclear derivative in :math:`H_{ep}` in the basis of these displacements. 

The matrix element describing the coupling between two nuclear wavefunctions which differ by replacing a single phonon :math:`\mathbf{q}j` with a phonon :math:`\mathbf{q}j'` therefore reads

.. math::
    
    g^{\mathbf{q}f,\mathbf{q}j} = -\frac{\hbar^2}{M}\sum_i \braket{\chi_{\mathbf{q}f}|\mathbf{e}_{\mathbf{q}j} \nabla_{\mathbf{R}}\chi_{\mathbf{q}i}} \cdot \sum_{s\mathbf{k}}w_{\mathbf{k}}\sum_{\nu<\epsilon_F}\sum_{\nu'>\epsilon_F} \braket{\psi_{s\mathbf{k+q}\nu'}|\mathbf{e}_{\mathbf{q}j}\nabla_{\mathbf{R}}\psi_{s\mathbf{k}\nu}}  ,

where Brillouin sampling weights :math:`w_{\mathbf{k}}` have been accounted for and final state indices follow from momentum conservation. Utilizing the ladder operator technique for harmonic oscillator eigenstates [Mahan2013]_, [Head-Gordon1992]_ we find

.. math::
   :label: app-HO-braket

   \braket{\chi_{\mathbf{q}f}| \mathbf{e}_{\mathbf{q}i}\nabla_{\mathbf{R}}\chi_{\mathbf{q}i}} &= \sqrt{\frac{\omega_{\mathbf{q}i} M}{2\hbar}} (\sqrt{n_{\mathbf{q},f}}\hat{\delta}_{\mathbf{q}f,\mathbf{q},i+1}- \sqrt{n_{\mathbf{q},f}+1}\hat{\delta}_{\mathbf{q}f,\mathbf{q},i-1}) .

As a result transitions can only occur between vibrational states with :math:`f=i\pm1`. If we assume that we start out in the first vibrationally excited state and are interested in decay to the ground state we arrive at following expression for the matrixelement

.. math::
   :label: app-matrixelement

    g^{\mathbf{q}j} &= \sqrt{\frac{\hbar^3\omega_{\mathbf{q}j}}{2M}}\sum_{s\mathbf{k}}w_{\mathbf{k}}  \sum_{\nu<\epsilon_F}\sum_{\nu'>\epsilon_F} g_{s\mathbf{k+q}\nu',s\mathbf{k}\nu}^{\mathbf{q}j} .
    
with

.. math::
   :label: matrixelement2

    g_{s\mathbf{k+q}\nu',s\mathbf{k}\nu}^{\mathbf{q}j} = \braket{\psi_{s\mathbf{k+q}\nu'}|\mathbf{e}_{\mathbf{q}j}\nabla_{\mathbf{R}}\psi_{s\mathbf{k}\nu}}.

Therefore, the excitation or deexcitation probability of a phonon :math:`\mathbf{q}j` is given by the coupling of all possible electronic excitations that do not violate momentum conservation. 

Often it is more desirable to express the matrix elements :math:`g_{s\mathbf{k+q}\nu',s\mathbf{k}\nu}^{\mathbf{q}j}` in a different form, namely in terms of nuclear derivatives of the unperturbed Hamiltonian and overlap matrices. In the case of a local atomic orbital representation of the wavefunctions, such a representation was derived by Head-Gordon and Tully [Head-Gordon1992]_ as well as Savrasov and Savrasov [Savrasov1996]_ in the current context. To review this derivation we express the wavefunctions in terms of a non-orthogonal, incomplete basis set

.. math::
  :label: app-lcao

  \ket{\psi_{s\mathbf{k}i}} = \sum_a c_{s\mathbf{k}i}^a \ket{\phi_{\mathbf{k}}^a}
     
such that following generalized eigenvalue matrix equation is satisfied (In the following intermediate results we will drop the spin index):

.. math::
  :label: app-eigenvalue

  \mathbf{H}_{\mathbf{k}} \mathbf{c}_{\mathbf{k}i} = \mathbf{S}_{\mathbf{k}} \mathbf{c}_{\mathbf{k}i}\epsilon_{\mathbf{k}i}
      
with

.. math::
    
    H_{\mathbf{k}}^{ab} = \braket{\phi_{\mathbf{k}}^b|\hat{H}_e|\phi_{\mathbf{k}}^a}
       
and

.. math::

    S_{\mathbf{k}}^{ab} = \braket{\phi_{\mathbf{k}}^b|\phi_{\mathbf{k}}^a} .

Substituting :eq:`app-lcao` into the Dirac-braket of eq. :eq:`app-matrixelement` gives

.. math::
   :label: app-braket

    \braket{\psi_{s\mathbf{k+q}\nu'}|\mathbf{e}_{\mathbf{q}j}\nabla_{\mathbf{R}}\psi_{s\mathbf{k}\nu}} = \sum_{a,b} c^b_{s\mathbf{k+q}\nu'}\braket{\phi_{\mathbf{k+q}}^b|\mathbf{e}_{\mathbf{q}j}\nabla_{\mathbf{R}}\phi_{\mathbf{k}}^a}c^a_{s\mathbf{k}\nu}  +c^b_{s\mathbf{k+q}\nu'}\braket{\phi_{\mathbf{k+q}}^b|\phi_{\mathbf{k}}^a}\mathbf{e}_{\mathbf{q}j}\nabla_{\mathbf{R}}c^a_{s\mathbf{k}\nu}

The second term of this equality can be reexpressed as follows using the nuclear derivative of eq. :eq:`app-eigenvalue`

.. math::
   :label: app-HandS-deriv

    c^b_{s\mathbf{k+q}\nu'}S_{\mathbf{k+q}\mathbf{k}}^{ab}\mathbf{e}_{\mathbf{q}j}\nabla_{\mathbf{R}} c^a_{s\mathbf{k}\nu} = \frac{1}{\epsilon_{\mathbf{k+q}\nu'}-\epsilon_{\mathbf{k}\nu}} c^b_{s\mathbf{k+q}\nu'}\left(H_{\mathbf{k+q},\mathbf{k}}^{\mathbf{q}j,ab}-\epsilon_{\mathbf{k}\nu} S_{\mathbf{k+q},\mathbf{k}}^{\mathbf{q}j,ab} \right) c^a_{s\mathbf{k}\nu}
         
with

.. math::

    S_{\mathbf{k+q},\mathbf{k}}^{ab} = \braket{\phi_{\mathbf{k+q}}^b|\phi_{\mathbf{k}}^a}

and the nuclear derivative of an operator in matrix representation defined as
          
.. math::

    O_{\mathbf{k+q},\mathbf{k}}^{\mathbf{q}j,ab} = \mathbf{e}_{\mathbf{q}j}\nabla_{\mathbf{R}}\braket{\phi_{\mathbf{k+q}}^b|O|\phi_{\mathbf{k}}^a} .

Inserting this result into eq. :eq:`app-matrixelement` and changing to matrix representation we find

.. math::
   :label: app-g-HandS

    g_{s\mathbf{k+q}\nu',s\mathbf{k}\nu}^{\mathbf{q}j} \approx\frac{ \mathbf{c}_{s\mathbf{k+q}\nu'}\left(\mathbf{H}_{\mathbf{k+q},\mathbf{k}}^{\mathbf{q}j}-\epsilon_{F} \mathbf{S}_{\mathbf{k+q},\mathbf{k}}^{\mathbf{q}j} \right) \mathbf{c}_{s\mathbf{k}\nu}} {\epsilon_{s\mathbf{k+q}\nu'}-\epsilon_{s\mathbf{k}\nu}}
            
In eq. :eq:`app-g-HandS` we assume that the energy difference between states :math:`\epsilon_{s\mathbf{k}\nu}` and :math:`\epsilon_{s\mathbf{k+q}\nu'}` is very small and that both are close to the Fermi level and can be replaced therewith [Head-Gordon1992]_. 

Unfortunately the above terms can not straightforwardly be extracted from an electronic structure calculation in reciprocal space representation due to the dependence on different wavevectors :math:`\mathbf{k}` and :math:`\mathbf{k+q}` and :math:`\braket{\mathbf{k}|\mathbf{H}|\mathbf{k}'}\neq \delta_{\mathbf{k}\mathbf{k}'}\cdot\braket{\mathbf{k}|\mathbf{H}|\mathbf{k}}`. In the case of a plane-wave representation a solution to this problem has been given by Trail *et al.* with coupling elements given in terms of eq. :eq:`app-matrixelement2` [Trail2001]_. In the current case, we can approach this by transforming the Hamiltonian and overlap matrices themselves to their real-space representations using a Fourier sum over wavevectors within the first Brillouin zone:

.. math::
   :label: app-Hreal

    \mathbf{H}(|\mathbf{N}-\mathbf{N}'|)=\sum_{\mathbf{k}} w_{\mathbf{k}} e^{i\mathbf{k}\cdot(\mathbf{R}_{\mathbf{N}'}-\mathbf{R}_{\mathbf{N}})}\mathbf{H}_{\mathbf{k}}

and
             
.. math::
   :label: app-Sreal

    \mathbf{S}(|\mathbf{N}-\mathbf{N}'|)= \sum_{\mathbf{k}} w_{\mathbf{k}} e^{i\mathbf{k}\cdot(\mathbf{R}_{\mathbf{N}'}-\mathbf{R}_{\mathbf{N}})} \mathbf{S}_{\mathbf{k}} .

In eq. :eq:`app-Hreal` and :eq:`app-Sreal`, :math:`\mathbf{N}` and :math:`\mathbf{N}'` are integer index-vectors that define the absolute position of a unit cell image in the crystal volume :math:`\Omega` and :math:`w_{\mathbf{k}}` is the normalization weight of wavevector :math:`\mathbf{k}`. We have furthermore defined crystal translation vectors as 

.. math::
    \mathbf{R}_\mathbf{N} = \mathbf{R}_0 + \mathbf{T}_\mathbf{N},

where :math:`\mathbf{R}_0` is the center of the principal unit cell and :math:`\mathbf{T}_\mathbf{N}` is a lattice translation vector to the unit cell image with index vector :math:`\mathbf{N}` [Soler2002]_, [Blum2009]_.

From the above results :math:`\mathbf{H}^{\mathbf{q}j}(|\mathbf{N}-\mathbf{N}'|)` and :math:`\mathbf{S}^{\mathbf{q}j}(|\mathbf{N}-\mathbf{N}'|)` can be directly calculated using the finite-difference method or Density Functional Perturbation Theory (Coupled Perturbed Kohn-Sham) [Savrasov1996]_, [Baroni2001]_. With the definition of the inverse transformation

.. math::
    \mathbf{H}^{\mathbf{q}j}_{\mathbf{k},\mathbf{k}'}=\sum_{\mathbf{N},\mathbf{N}'\in\Omega} e^{i(\mathbf{k}'\cdot\mathbf{R}_{\mathbf{N}'}-\mathbf{k}\cdot\mathbf{R}_{\mathbf{N}})} \mathbf{H}^{\mathbf{q}j}(|\mathbf{N}-\mathbf{N}'|)

and
                
.. math::
    \mathbf{S}^{\mathbf{q}j}_{\mathbf{k},\mathbf{k}'}= \sum_{\mathbf{N},\mathbf{N}'\in\Omega} e^{i(\mathbf{k}'\cdot\mathbf{R}_{\mathbf{N}'}-\mathbf{k}\cdot\mathbf{R}_{\mathbf{N}})} \mathbf{S}^{\mathbf{q}j}(|\mathbf{N}-\mathbf{N}'|)

we may finally express the coupling element of eq. :eq:`app-g-HandS` as

.. math::                  
    & g_{s\mathbf{k+q}\nu',s\mathbf{k}\nu}^{\mathbf{q}j} \approx\sqrt{\frac{\hbar}{M}}\frac{1}{\epsilon_{\mathbf{k+q}\nu'}-\epsilon_{\mathbf{k}\nu}}\sum_{\mathbf{N,N}'}e^{i((\mathbf{k+q})\cdot\mathbf{R}_{\mathbf{N}'}-\mathbf{k}\cdot\mathbf{R}_{\mathbf{N}})} \cdot \\ \nonumber &  \cdot \mathbf{c}_{s\mathbf{k+q}\nu'}\left(\mathbf{H}^{\mathbf{q}j}(|\mathbf{N}-\mathbf{N}'|)-\epsilon_{F} \mathbf{S}^{\mathbf{q}j}(|\mathbf{N}-\mathbf{N}'|) \right) \mathbf{c}_{s\mathbf{k}\nu}


Derivation using first order time dependent perturbation theory
----------------------------------------------------------------

Our starting point is given by eqs. :eq:`app-HO` and :eq:`app-matrixelement`. 
The harmonic nuclear wavefunctions are perturbed by a potential :math:`H_{\mathrm{ep}}`, 
which we assume to be harmonic in time, namely with the frequency :math:`\omega_{\mathbf{q}j}` of the vibration itself. 
This can be understood by considering that at the equilibrium position of the 
nuclei the friction effects due to the electronic structure are zero, 
while they are strongest at the maximum displacement. 
The effect, as already discussed, is that electronic and vibrational states hybridize. 
Therefore all 'vibronic' states that can be excited by the phonon frequency :math:`\omega_{\mathbf{q}j}` 
will contribute to the decay of the phonon. 
As a result, in the following we will not distinguish between energies of electronic or vibrational states. 
The potential introduces mixing between the nuclear eigenfunctions :math:`\chi`  characterized with coefficients :math:`c_n(t)`:

.. math::

    \ket{\chi(t)} = \sum_n c_n(t) \operatorname{e}^{-iE_n t/\hbar} \ket{\chi_n}

Expanding the coefficients and wavefunctions in a perturbation series and truncating after the first perturbation term, we can define the probability of excitation between two states (both differing by occupation change between phonons :math:`\mathbf{q}j` to :math:`\mathbf{q}f`) by the absolute square of the first order coefficients:

.. math::

    &P_{\mathbf{q}j\rightarrow \mathbf{q}f}(t)=|c^{(1)}_{jf}(t)|^2 = \\ \nonumber & =\frac{1}{\hbar^2}|\braket{\chi_{\mathbf{q}f}|H_{\mathrm{ep}}|\chi_{\mathbf{q}j}}|^2 \left( \frac{\sin((\Delta\omega-\omega_{\mathbf{q}j})t/2)}{(\Delta\omega-\omega_{\mathbf{q}j})} \right)^2 ,

with :math:`\Delta\omega` as energy difference between all vibronic states generated 
by the perturbation. Due to our previous assumption of harmonic oscillator states every 
phonon state can only decay into the next lowest state with the result that the overall 
probability of finding state :math:`\mathbf{q}j` being the excitation probability into this state from the next lowest state

.. math::
   :label: app-probab

   P_{\mathbf{q}j}(t) =|\braket{\chi_{\mathbf{q}f}|H_{\mathrm{ep}}|\chi_{\mathbf{q}i}}|^2 \left( \frac{\sin((\Delta\epsilon-\hbar\omega_{\mathbf{q}j})t/2\hbar)}{(\Delta\epsilon-\hbar\omega_{\mathbf{q}j})} \right)^2

At this point we can introduce our previous result for the coupling matrix element by inserting eq. :eq:`app-matrixelement` into eq. :eq:`app-probab`

.. math::
   :label: app-probab2
   
   P_{\mathbf{q}j}(t) = \frac{\hbar^3\omega_{\mathbf{q}j}}{2M}\sum_{s\mathbf{k}}w_{\mathbf{k}} \sum_{\nu<\epsilon_F}\sum_{\nu'>\epsilon_F}|g_{s\mathbf{k+q}\nu',s\mathbf{k}\nu}^{\mathbf{q}j}|^2 \left( \frac{\sin((\Delta\epsilon-\hbar\omega_{\mathbf{q}j})t/2\hbar)}{(\Delta\epsilon-\hbar\omega_{\mathbf{q}j})} \right)^2

We are interested in the decay rate in the limit of small perturbing frequencies or, equally, 
the limit of long time (:math:`t>>1/(\Delta\omega-\omega_{\mathbf{q}j})`). 
Taking this limit we can assume exponential decay and we arrive at the well known Fermi's Golden Rule expression:

.. math::    
    
    \Gamma_{\mathbf{q}j} = \lim_{t\rightarrow\infty} \frac{P_{\mathbf{q}j}(t)}{t} = \frac{\pi\hbar^2\omega_{\mathbf{q}j}}{M}\sum_{s\mathbf{k}}w_{\mathbf{k}} \sum_{\nu<\epsilon_F}\sum_{\nu'>\epsilon_F}|g_{s\mathbf{k+q}\nu',s\mathbf{k}\nu}^{\mathbf{q}j}|^2 
      \delta(\Delta\epsilon-\hbar\omega_{\mathbf{q}j})

However, even without taking the limit of :math:`t\rightarrow\infty`, applying the :math:`\sin`-function 
in eq. :eq:`app-probab2` with a reasonably chosen time would yield a similarly narrow 
integration kernel for the typical range of vibrational frequencies found for atomic or molecular adsorbates. 
Finally, we can generalize to non-integer finite-temperature occupations and arrive at:

.. math::
    :label: final-gamma
    
    \Gamma_{\mathbf{q}j} = \frac{\pi\hbar^2\omega_{\mathbf{q}j}}{M}\sum_{s\mathbf{k}}\sum_{\nu}\sum_{\nu'} (f(\epsilon_{s\mathbf{k+q}\nu'})-f(\epsilon_{s\mathbf{k}\nu}))|g_{s\mathbf{k+q}\nu',s\mathbf{k}\nu}^{\mathbf{q}j}|^2 
          \delta(\epsilon_{s\mathbf{k+q}\nu'}-\epsilon_{s\mathbf{k}\nu}-\hbar\omega_{\mathbf{q}j}) .


With :math:`\hbar\omega_{\mathbf{q}j}` being the perturbing energy for which we apply the delta function in :eq:`final-gamma` we can re-express the above equation as

.. math::
    
    \Gamma_{\mathbf{q}j} = \frac{\pi\hbar}{M}\sum_{s\mathbf{k}}\sum_{\nu}\sum_{\nu'} (f(\epsilon_{s\mathbf{k+q}\nu'})-f(\epsilon_{s\mathbf{k}\nu}))|g_{s\mathbf{k+q}\nu',s\mathbf{k}\nu}^{\mathbf{q}j}|^2 \cdot (\epsilon_{s\mathbf{k+q}\nu'}-\epsilon_{s\mathbf{k}\nu}) \cdot\delta(\epsilon_{s\mathbf{k+q}\nu'}-\epsilon_{s\mathbf{k}\nu}-\hbar\omega_{\mathbf{q}j}) .

For details on how this equation is used in this software, see :ref:`theory`.


References
----------

.. [Maurer2015] R.J. Maurer, M. Askerka, and J. Tully, *in preparation*  (2015)
.. [Mahan2013] G.D. Mahan, *Many-particle physics*, Springer Science & Business Media (2013) 
.. [Savrasov1996] `S.` Savrasov, and D. Savrasov, *Phys. Rev. B* **54**, 16487-16501 (1996) 
.. [Baroni2001] `S.` Baroni, S. de Gironcoli, and A. Dal Corso, *Rev. Mod. Phys.* **73**, 515-562 (2001)
.. [Blum2009] `V.` Blum, R. Gehrke, F. Hanke, P. Havu, V. Havu, X. Ren, K. Reuter, adM. Scheffler, *Comp. Phys. Comm.* **180**, 2175-2196 (2009)
.. [Soler2002] J.M. Soler, E. Artacho, J.D Gale, A. Garcia, J. Junquera, P. Ordejon, and D. Sanchez-Portal, *J. Phys.: Condens. Matter* **14**, 2745 (2002)
.. [Trail2001] `J.` Trail, M. Graham, and D. Bird, *Comp. Phys. Comm.* **137**, 163-173 (2001)
