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
from typing import Any, Union

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
def dist(a, e):
    return a * (1 + (e ** 2) / 2)


# orbital average stellar flux [1]
# a = semi-major axis (AU), e = eccentricity, L = luminosity (solar units)
def flux(a, e, L):
    return L / (a ** 2 * math.sqrt(1 - e ** 2))


# orbital average equilibrium temperature [1]
# a = semi-major axis (AU), e = eccentricity, A = albedo, L = luminosity (solar units)
def teq(a, e, A, L):
    return 278.5 * ((1 - A) * L / a ** 2) ** (1 / 4.) * (2 * math.sqrt(1 + e) / math.pi) * \
           special.ellipe(2 * e / (1 + e))


# effective thermal distance [1]
# a = semi-major axis (AU), e = eccentricity
def reff(a, e):
    return a * ((2 * math.sqrt(1 + e) / math.pi) * scipy.special.ellipe(2 * e / (1 + e))) ** (-2.)


def edist(a, e, Eo):
    return a * ((((e ** 2) + 2) * Eo + e * (e * math.cos(Eo - 4)) * (
        math.sin(Eo))) / (2 * (Eo - e * math.sin(Eo))))


def eflux(a, e, L, Eo):
    return (L / a ** 2) * ((2 * math.atan(math.sqrt(1 + e / 1 - e)) * (math.tan(1 / 2) * Eo)) / (
            (Eo - math.sin(e) * Eo) * math.sqrt(1 - e ** 2)))


def eteq(a, e, A, L, Eo):
    return 278.5 * ((L * (1 - A) / a ** 2) ** (1 / 4)) * ((2 * 1j * math.sqrt(1 - e)) / (Eo - e * math.sin(Eo))) * (
            special.ellipkinc((1 + e) / (1 - e), 1j * math.asinh(math.tan(1 / 2) * Eo)) - special.ellipeinc(
        (1 + e) / (1 - e), 1j * math.asinh(math.tan(1 / 2) * Eo)) - 1j * math.sqrt((1 - e * math.cos(Eo)) / (1 - e)) *
            math.tan(1 / 2) * Eo)


def preff(Do, rp):
    return rp((3 * Do ** 4 + 10 * Do ** 2 + 15) / (5 * (Do ** 2 + 3)))


def pflux(L, rp):
    return (L / rp ** 2) * ((3 * math.atan(Do)) / ((Do ** 2 + 3) * Do))


def pteq(A, L, Do):
    return 278.5 * (((L * (1 - A)) / (rp ** 2)) ** (1 / 4)) * (
            (3 * (Do * math.sqrt(Do ** 2 + 1) + math.asinh(Do))) / (2 * (Do ** 2 + 3) * Do))


def hreff(e, Ho):
    return rp * ((e * (e * math.cosh(Ho) - 4) * math.sinh(Ho) + ((e ** 2) + 2) * Ho) / (
            2 * (e * math.sinh(Ho) - Ho) * (e - 1)))


def hflux(e, L, Ho, rp):
    return (L / rp ** 2) * ((2 * ((e - 1) ** 2) * (((e ** 2) - 1) ** (1 / 2)) / (e * math.sinh(Ho) - Ho)) * math.atan(
        ((e + 1) * math.tanh(1 / 2) * Ho) / (math.sqrt((e ** 2) - 1))))


def hteq(e, A, L, Ho):
    return 278.5 * ((L * (1 - A)) ** (1 / 4)) * (
            ((-2 * math.sqrt(-1) * (e - 1)) / (e * math.sinh(Ho) - Ho)) * special.ellipkinc(
        (1 / 2) * math.sqrt(-1) * Ho, ((2 * e) / (e - 1))))


def lreff(rp, Lo):
    return rp * ((math.sinh(2 * Lo) + 2 * Lo) / (4 * math.sinh(Lo)))


def lflux(L, rp, Lo):
    return (L / rp ** 2) * (2 * math.atan(math.tanh((1 / 2) * Lo)) / math.sinh(Lo))


def lteq(A, L, Lo, rp):
    return 278.5 * ((L * (1 - A) / rp ** 2) ** (1 / 4)) * (
            -2j / math.sinh(Lo) * special.ellipkinc((1 / 2) * 1j * Lo, 2))
