# 22.312, PSet10
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 13-7: Calculation of CPR for a BWR hot channel

from iapws import IAPWS97 as Steam
from scipy.optimize import fsolve
from scipy.integrate import quad as integral
from pylab import *

NPOINTS = 50
KELVIN = 273.15
# Given constants
Q1REF = 104.75              # kW/m
ALPHA = 1.96                # <dimensionless>
G = 1569.5                  # kg/s/m^2
KGSM2_TO_LBHFT2 = 737.3E-6  # kg/s/m^2 -> Mlbm/h/ft^2
G_WEIRD = G*KGSM2_TO_LBHFT2 # Mlbm/h/ft^2
L = 3.588                   # m
D = 0.0112                  # m
AF = 9.718E-3               # m^2
ACH = 1.42E-4               # m^2
NCH = 74                    # channels
P = 7.14                    # MPa
TIN = 278.3 + KELVIN        # K
# Hench-Gillis correlation
J = 1.032
FP = -1.66E-3
A = 0.5*G_WEIRD**-0.43
B = 165 + 115*G_WEIRD**2.3
# Steam tables
water_in = Steam(P=P, T=TIN)
hin = water_in.h
sat_water = Steam(P=P, x=0)
hf = sat_water.h
sat_vapor = Steam(P=P, x=1)
hfg = sat_vapor.h - hf

def hench_gillis(lboil):
	"""The Hench-Gillis correlation for critical quality
	
	Parameter:
	----------
	lboil:  float, m; boiling length of the BWR channel
	
	Returns:
	--------
	"""
	zed = pi*D*NCH*lboil/AF
	xc = A*zed/(B + zed) * (2 - J) + FP
	return xc

def q1(z):
	"""Given linear power profile
	
	Parameter:
	----------
	z:      float, m; axial location. z=0 is the midplane
	
	Returns:
	--------
	q':     float, kW/m; linear power at that location
	"""
	return Q1REF*exp(-ALPHA*(z/L + 0.5))*cos(pi/L*z)

print("A = {:.2f}, \t B = {:.1f}".format(A, B))
mdot = G*ACH
print("mdot:     {:.3f} kg/s".format(mdot))
xein = (hin - hf)/hfg
print("Xe,in:    {:.3f}".format(xein))
# Numerically solve for the onset of nucleate boiling
f = lambda z: integral(q1, -L/2, z)[0]/(mdot*hfg) + xein
zonb = fsolve(f, 0.0)[0]
print("zonb:     {:.3f} m ({:.2f} m from inlet)".format(zonb, L/2 + zonb))

# Finally, find the critical power ratio
# Critical power is the point where the Hench-Gillis correlation for X_cr meets
# the X_e curve. CPR is the factor times q1(z) necessary to make that occur.
zrange = linspace(zonb, L/2, NPOINTS)
xcr_vals = array([hench_gillis(z - zonb) for z in zrange])

def _xe(cpr, z):
	_q1 = lambda _z: cpr*q1(_z)
	qdot = integral(_q1, zonb, z)[0]
	h = hf + 1.5E3*qdot/(mdot*hfg)
	#print(h, h-hf, hf)
	return (h - hf)/hfg

def check_cpr(cpr, tol = 1E-3):
	"""Multiply `cpr` times q1(z) and evaluate whether the
	curves intersect within tolerance.
	
	Parameter:
	----------
	cpr:        float; Critical Power Ratio.
	tol:        float; absolute tolerance
	
	Returns:
	--------
	Boolean; whether the curves intersect
	"""
	xe_vals = array([_xe(cpr, z) for z in zrange])
	plot(zrange, xe_vals, label = "$X_e$ ({:.0%})".format(cpr))
	diffs = xcr_vals - xe_vals
	if min(diffs[1:]) < tol:
		return True
	else:
		return False

def get_cpr(step = 0.03):
	guess = 1.0
	while True:
		if check_cpr(guess):
			return guess
		else:
			guess += step

cpr = get_cpr()
txt = "CPR = {:.0%}".format(cpr)
print(" -> " + txt)

plot(zrange, xcr_vals, "k-", label = "$X_c$ (Hench-Gillis)")
xe_vals = array([_xe(1, z) for z in zrange])
xlim([-L/2, +L/2])
ylim([0, 1.0/3])
title(txt)
legend()
grid()
show()
