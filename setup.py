import os
from setuptools import setup

# Utility function to read the README file.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "vibcooling",
    version = "0.0.1",
    author = "Mikhail Askerka & Reinhard J. Maurer",
    author_email = "mikhail.askerka@yale.edu, reinhard.maurer@yale.edu",
    description = ("
            This package contains routines and scripts to calculate 
            electron-phonon coupling, friction tensors, and vibrational lifetimes 
            from different quantum chemistry codes.
            "),
    license = "TBA",
    keywords = "vibrational cooling, quantum chemistry",
    url = "TBA",
    package_dir = {'': 'src'},
    py_modules=[
        'parser', 
        ],
    long_description=read('README'),
)
