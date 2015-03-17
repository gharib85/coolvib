import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext 

# Utility function to read the README file.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

package_name = 'vibcooling'

__version__ = '0.0.1'

__all__ = [
        'parser',
        'tools',
        ]

packages = [ package_name ]
for package in __all__:
    packages.append(package_name + '.' + package)

ext_modules = [
           Extension("vibcooling.parser.siesta_mod",["src/parser/siesta_mod.pyx"]) 

        ]

setup(
    name = package_name,
    version = __version__,
    author = "Mikhail Askerka & Reinhard J. Maurer",
    author_email = "mikhail.askerka@yale.edu, reinhard.maurer@yale.edu",
    description = ("This package contains routines and scripts to calculate \
            electron-phonon coupling, friction tensors, and vibrational lifetimes \
            from different quantum chemistry codes. \
            "),
    license = "TBA",
    keywords = "vibrational cooling, quantum chemistry",
    url = "TBA",
    cmdclass = {'build_ext' : build_ext},
    package_dir = {package_name: 'src'},
    packages = packages,
    ext_modules= ext_modules,
    long_description=read('README'),
)
