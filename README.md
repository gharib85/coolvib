# README #

### What does coolvib do? ###

Coolvib uses Fermi's Golden rule to calculate vibrational relaxation rates 
and lifetimes due to nonadiabatic effects on metal surfaces. The code uses 
precomputed electronic structure information from software packages such as 
FHI-aims and SIESTA.

For more details, have a look at the documentation at https://maurergroup.github.io/coolvib/

### Dependencies ###

* Python 2.7.x
* Numpy >=1.6
* Scipy >=0.12
* matplotlib>=1.2
* Atomic Simulation Environment [ ASE ](https://wiki.fysik.dtu.dk/ase/)

### Installation ###

* After installing all dependencies, just issue make 
and include build/ into PYTHONPATH

* export PYTHONPATH=<path to coolvib>/coolvib/build/lib.<arch>-x86_64-<py-version>:$PYTHONPATH
* export PATH=<path to coolvib>/coolvib/build/lib.<arch>-x86_64-<py-version>/coolvib/tools:$PATH

## Licensing ##

coolvib is licensed under the GNU General Public License, version 3 (gnu.org/licenses/gpl.html)

### Who do I talk to? ###

r.maurer@warwick.ac.uk

### Contributors ###
Reinhard Maurer,  
Mikhail Askerka
