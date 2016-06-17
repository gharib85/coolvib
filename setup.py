#    This file is part of coolvib
#
#        coolvib is free software: you can redistribute it and/or modify
#        it under the terms of the GNU General Public License as published by
#        the Free Software Foundation, either version 3 of the License, or
#        (at your option) any later version.
#
#        coolvib is distributed in the hope that it will be useful,
#        but WITHOUT ANY WARRANTY; without even the implied warranty of
#        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#        GNU General Public License for more details.
#
#        You should have received a copy of the GNU General Public License
#        along with coolvib.  If not, see <http://www.gnu.org/licenses/>.

import os
#from distutils.core import setup
from distutils.core import setup, Command
from distutils.extension import Extension
from Cython.Distutils import build_ext 
import numpy as np

# Utility function to read the README file.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

package_name = 'coolvib'

__version__ = '0.1'

__all__ = [
    'parser',
    'tools',
    'routines',
    'convenience',
    'test',
    ]

packages = [ package_name ]
for package in __all__:
    packages.append(package_name + '.' + package)

class test(Command):
    def __init__(self, dist):
        Command.__init__(self, dist)
        self.sub_commands = ['build']

ext_modules = [
           Extension("coolvib.parser.siesta_mod",["coolvib/parser/siesta_mod.pyx"],
               include_dirs=[np.get_include()]), 
           Extension("coolvib.routines.build_G",["coolvib/routines/build_G.pyx"],
               include_dirs=[np.get_include()]) 
        ]

setup(
    name = package_name,
    version = __version__,
    url = "www.damaurer.at",
    author = "Reinhard J. Maurer",
    author_email = "reinhard.maurer@yale.edu",
    description = ("This package contains routines and scripts to calculate \
            electron-phonon coupling, friction tensors, and vibrational lifetimes \
            from different quantum chemistry codes. \
            "),
    license = "GNU General Public License",
    keywords = "vibrational cooling, quantum chemistry",
    cmdclass = {'build_ext' : build_ext,
        'test' : test},
    packages = packages,
    package_dir = {package_name: package_name},
    platforms='linux',
    ext_modules= ext_modules,
    long_description=read('README.md'),
)
