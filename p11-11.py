# 22.312, PSet08
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 11-11: HEM Pressure loss

import models
import two_phase
from math import pi, sqrt
from scipy.integrate import quad as integral
from iapws import IAPWS97 as Steam

# Given constants
L = 3                   # m
AFLOW = 1.5E-4          # m^2
DH = 2*sqrt(AFLOW/pi)   # m
MDOT = 0.29             # kg/s
G = MDOT/AFLOW          # kg/s/m^2
P = 7.2                 # MPa
X1 = 0.15
# Steam tables
water = Steam(P=P, x=0)
rhof = water.rho
mu = water.mu
vapor = Steam(P=P, x=1)
rhog = vapor.rho

print("Steam tables: rhof = {:.1f} kg/m^3, {:.1f} kg/m^3, \
mu = {:.2e} Pa-s".format(rhof, rhog, mu))
print("Mass flux:    {:.0f} kg/s/m^2".format(G))
print("\nPart 1) Adiabatic channel with Xin = {}".format(X1))
alpha1 = two_phase.alpha(X1, rhof, rhog)
rhom1 = rhog + (1 - alpha1)*(rhof - rhog)
print("\tMixture:      alpha = {:.3f}, rhom = {:.1f} kg/m^3".format(alpha1, rhom1))
re = G*DH/mu
print("\tReynolds:     {:.2e}".format(re))
fff = models.mcadams(re)
print("\tFric. fac.:   {:.3f}".format(fff))
dp1 = rhom1*9.81*L + fff*L/DH*G**2/(2*rhom1)
print("\t -> DP = {:.3} kPa".format(-dp1/1000))

print("\nPart 2) Uniform heat flux")
def x(z, x0 = 0):
	"""Quality of the mixture on [0, L]
	Requires saturated liquid entering the tube.
	
	Parameters:
	-----------
	z:          float, m; distance from tube inlet at bottom
	x0:         float; inlet quality
				[Default: 0]
	"""
	return x0 + z/L*X1
	
def alpha(z):
	"""Void fraction of the mixture on [0, L]
	
	Parameter:
	----------
	z:      float, m; distance from the inlet at the bottom
	"""
	return two_phase.alpha(x(z), rhof, rhog)

def rho(z):
	"""Mixture density on [0, L]
	
	Parameter:
	----------
	z:      float, m; distance from the inlet at the bottom
	"""
	return two_phase.mixture_density(rhof, rhog, alpha(z))

grav = lambda z: rho(z)*9.81
dp_grav = integral(grav, 0, L)[0]
print("\tDeltaP_grav: {:.2f} kPa".format(dp_grav/1000))
fric = lambda z: fff/DH*G**2/(2*rho(z))
dp_fric = integral(fric, 0, L)[0]
print("\tDeltaP_fric: {:.2f} kPa".format(dp_fric/1000))
dp_accl = G**2 * (1/rhom1 - 1/rhof)
print("\tDeltaP_accl: {:.2f} kPa".format(dp_accl/1000))
dp_total = dp_grav + dp_fric + dp_accl
print("\t -> DeltaP: {:.2f} kPa".format(-dp_total/1000))
