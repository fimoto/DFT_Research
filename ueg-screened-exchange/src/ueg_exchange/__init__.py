"""Screened second-order exchange in the uniform electron gas.

The JPA paper develops the spectral representation, diagram weights and
density limits. The PRB paper derives the spin-dependent low-density expansion.
See kernels_review/equation_map.json and spin/equation_map.json for equation labels.

Modules
-------
constants     units, normalization, asymptotic coefficients and spin moments
lindhard      imaginary-axis Lindhard response and RPA screened interactions
lorentzian    constant-width Lorentzian kernel and radial-energy asymptotes
dispersion    screening-pole approximation and leading low-density coefficient
highdensity   high-density logarithms and exchange self-energy partner
lowdensity    three-term low-density coefficients and spin dependence
ringsum       complementary RPA ring-energy checks
brm           direct BRM integral evaluator
"""
from . import constants, lindhard, lorentzian, dispersion, highdensity, lowdensity, ringsum
__version__ = "0.2.0"
