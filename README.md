The Planetary Habitability Library

	This Python module contains a collection of functions related to planetary
	habitability. This is a work in progress and the first version is planned
	for December 2019.

Functions (references in brackets):

	teq(a,e,A,L) - orbital average equilibrium temperature (K) [1]
	dist(a,e) - orbital average distance (AU) [1]
	flux(a,e,L) - orbital average stellar flux (solar units) [1]
	reff(a,e) - orbital effective thermal distance (AU) [1]

Definitions:
	
	a - orbital semi-major axis (AU)
	e - orbital eccentricity
	A - planetary bond albedo
	L - stellar luminosity (solar units)

Usage:

	Copy phl.py to your Python working directory. 
	
	Example:
	
	In : import phl
	In : phl.teq(1.0, 0.0, 0.3, 1.0)
	Out: 254.74150455519143
	
Requirements:

	Requires the Math and SciPy Python libraries.

References:

	[1] Méndez, Abel, and Edgard G. Rivera-Valentín. “The Equilibrium Temperature of
	Planets in Elliptical Orbits.” The Astrophysical Journal Letters 837, no. 1 (2017):
	L1. https://doi.org/10.3847/2041-8213/aa5f13.