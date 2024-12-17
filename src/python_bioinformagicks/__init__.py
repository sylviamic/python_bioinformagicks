'''
Top-level package for python bioinformagicks.
'''

__author__ = """Sylvia N. Michki"""
__email__ = 'sylvia.michki@gmail.com'
__version__ = '0.1.0'

from . import tools as tl
from . import plotting as pl

__all__ = [
    "tl",
    "pl",
]