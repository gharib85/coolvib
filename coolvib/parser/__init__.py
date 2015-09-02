"""
parser module

This module contains all necessary routines to read
precalculated data from quantum chemistry packages. 

There are routines that bundle the IO for certain purposes, such as 
:py:func:`coolvib.parser.aims.parse_aims_tensor`, but one can also 
use the individual routines.

"""

codes = {
        'aims': 'parse_aims',
        'siesta': 'parse_siesta',
        } 

code_type = {
        'aims' : 'local',
        'siesta' : 'local',
        }


from coolvib.parser.siesta import parse_siesta_tensor, parse_siesta_mode
from coolvib.parser.aims import parse_aims_tensor, parse_aims_mode
