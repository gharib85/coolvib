import os
from setuptools import setup

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
    package_dir = {package_name: 'src'},
    packages = packages, 
    long_description=read('README'),
)
