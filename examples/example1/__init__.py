"""
Example 1: Calculating a friction tensor from FHI-aims data: CO on Cu(100)
=========================================================================================

In this example we calculate the friction tensor from data previously 
calculated with `FHI-Aims <https://aimsclub.fhi-berlin.mpg.de/>`_.
We do this for CO adsorbed at a top site of a Cu(100) surface. The corresponding 
friction tensor is a (6x6) matrix.

In order to run this example you need to download the FHI-Aims dataset 
from following URL and extract it in the example folder

`Download Dataset <http://www.damaurer.at/software/coolvib/example1.tar.gz>`_


Code
-----

.. literalinclude:: ../examples/example1/example1.py

Detailed Explanation
---------------------

System Setup
^^^^^^^^^^^^

We start by setting up the system using an `ASE atoms object <https://wiki.fysik.dtu.dk/ase/ase/atoms.html#module-ase.atoms>`_
and we initialize our workflow_tensor class with

.. code-block:: python

    import coolvib
    model = coolvib.workflow_tensor(system, code='aims', active_atoms=[3,4])

, where *active_atoms* specifies the indices of the atoms we want to include in 
the calculation of the friction tensor. *Model* is an instance of the :py:class:`coolvib.convenience.workflow_tensor.workflow_tensor` class that guides us through the calculation and bundles all data and functions.

Reading input data
^^^^^^^^^^^^^^^^^^

After defining all keywords we can read the FHI-Aims data by 
using :py:func:`coolvib.convenience.workflow_tensor.workflow_tensor.read_input_data`.

The parser assumes that the input data is given as a set of directories with *eq* refering 
to the calculation at equilibrium position and *a3c0+* or *a3c0-* refering to displacements 
in positive and negative directions for atom3 and the cartesian x component (c0 = x, c1 = y, c2 = z).

Calculating the spectral function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is an optional step when calculating the friction tensor and at the moment only serves analysis purposes.

Using :py:func:`coolvib.convenience.workflow_tensor.workflow_tensor.calculate_spectral_function` we can calculate the spectral functions along all cartesian components and all combinations of cartesian components. This includes loops over all possible excitations and construction of the electron-phonon spectral function on a discrete grid. This discretization is determined by following keywords:

    *discretization_type*
        Specifies which broadening function is used to discretize the individual excitations, for more details see :py:func:`coolvib.routines.discretize_peak`

    *discretization_broadening*
        Specifies the width that is used to generate a finite width discrete representation of the peak and intensity (in eV)

    *discretization_length*
        Specifies the bin width of equidistant grid on which the spectral function is represented (in eV).

    *max_energy*
        Specifies the maximum extend beyond the fermi level for which the spectral function is calculated. The default is 3 eV. A safe value is 5 times the delta_function_width+perturbing_energy

The calculated spectral functions can be printed to file using :py:func:`coolvib.convenience.workflow_tensor.workflow_tensor.print_spectral_function`. This file can be post-processed and further analyzed using scripts from :py:mod:`coolvib.tools`.

.. image:: ./spectrum.png
   :align: center
   :width: 85%

Calculating the friction tensor and analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to finally arrive at the friction tensor one either has to integrate the spectral functions 
over a delta function (see :py:func:`coolvib.routines.evaluate_delta_function`). This is where the 
main approximations of the underlying method are applied. The width, position and shape of the 
delta function are controlled by following keywords:

    *delta_function_type*
        Specifies the type of delta function representation to be used

    *delta_function_width*
        Specifies the width of the delta function representation (in eV)

    *perturbing_energy*
        Specifies the position of the delta function representation. Typically this is chosen to be at the Fermi level under the quasi-static approximation (see :ref:`theory`)

Convergence with respect to the decay rates and lifetimes has to be validated for the above parameters.

The function :py:func:`coolvib.convenience.workflow_tensor.workflow_tensor.calculate_friction_tensor_from_spectrum` performs this integration and sets up the tensor. 

One can directly calculate the friction tensor by using :py:func:`coolvib.convenience.workflow_tensor.calculate_friction_tensor`, as it is done in this example

Analysis and printing of the tensor, its eigenvalues and eigenfunctions as well as the inverse of the eigenvalues (the principal lifetimes) can be done using :py:func:`coolvib.convenience.workflow_tensor.workflow_tensor.analyse_friction_tensor`.

Results
^^^^^^^

The expected results are 

Friction Tensor
      
             0                     1                     2                     3                     4                     5  
     0 3.1366e-02 7.7414e-07 7.7374e-06 -1.6675e-02 -2.7957e-07 -3.0235e-06 
     1 7.7414e-07 3.1365e-02 7.0913e-06 1.3343e-08 -1.6675e-02 -2.6273e-06 
     2 7.7374e-06 7.0913e-06 9.0953e-03 -2.1645e-06 -2.1907e-06 -2.8059e-03 
     3 -1.6675e-02 1.3343e-08 -2.1645e-06 9.6877e-03 -1.2865e-07 6.4511e-07 
     4 -2.7957e-07 -1.6675e-02 -2.1907e-06 -1.2865e-07 9.6872e-03 4.9371e-07 
     5 -3.0235e-06 -2.6273e-06 -2.8059e-03 6.4511e-07 4.9371e-07 1.6467e-03 
      
 
Friction Eigenvalues in 1/ps
         0          1          2          3          4          5
    0.0006     0.0006     0.0007     0.0100     0.0404     0.0404
Friction Eigenvectors
      
             0                     1                     2                     3                     4                     5  
     0 -1.3381e-01 4.5781e-01 -1.7908e-03 1.4193e-04 4.9189e-01 -7.2838e-01 
     1 -4.5780e-01 -1.3383e-01 -2.1490e-03 1.3969e-04 -7.2838e-01 -4.9190e-01 
     2 -1.4942e-03 6.0177e-04 3.1725e-01 -9.4834e-01 -5.3986e-05 -3.4853e-04 
     3 -2.4659e-01 8.4361e-01 -3.3172e-03 -3.1598e-04 -2.6695e-01 3.9527e-01 
     4 -8.4361e-01 -2.4659e-01 -3.9656e-03 -2.7458e-04 3.9528e-01 2.6694e-01 
     5 -5.1806e-03 2.2797e-03 9.4832e-01 3.1726e-01 1.5499e-05 1.2534e-04 
      
Principal lifetimes in ps
         0          1          2          3          4          5
 1566.7157  1566.3840  1412.3889    99.6617    24.7440    24.7431



"""
