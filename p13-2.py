# 22.312, PSet09
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 13-2: Factors affecting incipient superheat in a flowing system

import models
from math import pi, sqrt
from scipy.optimize import fsolve
from iapws import IAPWS97 as Steam

KELVIN = 273.15
# Given constants
P1 = 0.1013         # MPa
T2 = 290 + KELVIN   # K
HTC = 10000         # W/m^2-K
D = 0.02            # m
A = pi/4*D**2       # m^2
# Steam tables
water1 = Steam(P=P1, x=0)
vapor1 = Steam(P=P1, x=1)
water2 = Steam(T=T2, x=0)
vapor2 = Steam(T=T2, x=1)

def onb(water, vapor, htc, dt):
	"""Equation to solve for the temperature difference in order to find
	the point at which the onset of nucleate boiling occurs.
	
	Parameters
	-----------
	water:          IAPWS97 for saturated liquid
	vapor:          IAPWS97 for saturated vapor
	htc:            float, W/m^2-K; heat transfer coefficient
	dt:             float; dummy variable for the superheat
	
	Returns:
	--------
	hopefully 0
	"""
	hfg = 1000*(vapor.h - water.h)
	tsat = water.T
	return 2*sqrt(2*water.sigma*tsat*vapor.v*htc*dt/
	              (hfg*water.k)) - dt


print("\nPart 1: Atmospheric pressure")
onb1 = lambda dt: onb(water1, vapor1, HTC, dt)
dt1 = fsolve(onb1, 1)[0]
q21 = HTC*dt1
print("\t -> DeltaT:        {:.2f} K".format(dt1))
print("\t -> q'':           {:.2f} kW/m^2".format(q21/1000))

print("\nPart 2: 290 degreese Celsius")
onb2 = lambda dt: onb(water2, vapor2, HTC, dt)
dt2 = fsolve(onb2, 1)[0]
q22 = HTC*dt2
print("\t -> DeltaT:        {:.3f} K".format(dt2))
print("\t -> q'':           {:.1f} W/m^2".format(q22))

print("\nPart 3: Flow rate doubled from Part 2")
nu2 = HTC*D/water2.k
print("\tNu2 =              {:.1f}".format(nu2))
pr = water2.Prandt
print("\tPr:                {:.2f}".format(pr))
g2 = (nu2/(0.023*pr**0.4))**1.25 * water2.mu/D
print("\tG2:                {:.0f} kg/s/m^2".format(g2))
g3 = g2*2
print("\tG:                 {:.0f} kg/s/m^2".format(g3))
re3 = g3*D/water2.mu
print("\tRe:                {:.2e}".format(re3))
nu = models.dittus_boelter(re3, pr)
print("\tNu:                {:.1f}".format(nu))
h3 = nu*water2.k/D
print("\th:                 {:.2f} kW/m^2-K".format(h3/1000))
onb3 = lambda dt: onb(water2, vapor2, h3, dt)
dt3 = fsolve(onb3, 1)[0]
q23 = h3*dt3
print("\t -> DeltaT:        {:.3f} K".format(dt3))
print("\t -> q'':           {:.1f} W/m^2".format(q23))
