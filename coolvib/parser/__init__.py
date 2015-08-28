"""
parser module
"""

codes = {
        'aims': 'parse_aims',
        'siesta': 'parse_siesta',
        } 

code_type = {
        'aims' : 'local',
        'siesta' : 'local',
        }


from coolvib.parser.siesta import parse_siesta
from coolvib.parser.aims import parse_aims
