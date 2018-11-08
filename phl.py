#!/usr/bin/env python

"""
The Planetary Habitability Library

This Python module contains a collection of functions related to planetary habitability.
This is a work in progress and the first version is planned for December 2018.

Example:
	In : import phl
	In : phl.teq(1.0,0.0,0.3,1.0)
	Out: 254.74150455519143

References:

[1]	Méndez, Abel, and Edgard G. Rivera-Valentín. “The Equilibrium Temperature of
	Planets in Elliptical Orbits.” The Astrophysical Journal Letters 837, no. 1 (2017):
	L1. https://doi.org/10.3847/2041-8213/aa5f13.

"""

__author__ = "Abel Mendez"
__copyright__ = "PHL (2018)"
__credits__ = "PHL @ UPR Arecibo"
__license__ = "GPL"
__version__ = "0.0"
__maintainer__ = "Abel Mendez"
__email__ = "abel.mendez@upr.edu"
__status__ = "Production"

import math
from scipy import special

# orbital average distance [1]
# a = semi-major axis (AU), e = eccentricity
def dist(a,e):
	return a*(1+(e**2)/2)
	
# orbital average stellar flux [1]
# a = semi-major axis (AU), e = eccentricity, L = luminosity (solar units)
def flux(a,e,L):
	return L/(a**2*math.sqrt(1 - e**2))
	
# orbital average equilibrium temperature [1]
# a = semi-major axis (AU), e = eccentricity, A = albedo, L = luminosity (solar units)
def teq(a,e,A,L):
	To = 278.5
	return To*((1-A)*L/a**2)**(1/4.)*(2*math.sqrt(1+e)/math.pi)*\
	special.ellipe(2*e/(1+e))
	
# effective thermal distance [1]
# a = semi-major axis (AU), e = eccentricity
def reff(a,e):
	return a*((2*math.sqrt(1+e)/math.pi)*special.ellipe(2*e/(1+e)))**(-2.)