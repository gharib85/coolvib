"""
parser module
"""

codes = {
        'aims': 'parse_aims',
        'siesta': 'parse_siesta',
       ] 

import coolvib.parser.parser as parser

parse_aims = parser.parse_aims
parse_siesta = parser.parse_siesta
