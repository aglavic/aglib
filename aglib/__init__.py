'''
Module for stuff often used in scattering.
'''

from . import functions
from .functions import *
from .peakfinder import PeakFinder
from .mpfit import mpfit

__all__=['FF', 'res', 'mpfit', 'PeakFinder']+functions.__all__
