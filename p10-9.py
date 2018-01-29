# 22.312, PSet07
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 10.9
# AHTR

import models
from math import pi, log, sqrt
from scipy.optimize import fsolve

# Given parameters
DPFRIC = 200000         # Pa
TIN = 600               # deg C
TMAX = 1000             # deg C
L = 10                  # m; length of fuel
DFLAT = 0.03            # m; flat-to-flat hexagon with
DFLOW = 0.01            # m; coolant channel diameter
AFLOW = pi/4*DFLOW**2   # m^2; coolant flow area


class Material(object):
	"""Material with thermal hydraulic properties

	Parameters:
	-----------
	:param rho:     float, kg/m^3; density
	:param k:       float, W/m-K;  conductivity
	:param mu:      float, Pa-s;   viscosity
	:param c:       float, J/kg-K; heat capacity
	
	Attributes:
	-----------
	pr:             float; Prandtl number
	"""
	def __init__(self, rho, k, mu, c):
		self.rho = rho
		self.k = k
		self.mu = mu
		self.c = c
		if self.mu:
			self.pr = c*mu/k
		else:
			self.pr = None
	
	def get_reynolds(self, g, lc = DFLOW):
		"""Return the Reynolds number at some mass flux
		
		Parameters:
		-----------
		g:          float, kg/s/m^2; mass flux
		lc:         float, m; characteristic length of flow
					[Default: DFLOW from constants]
		"""
		return g*lc/self.mu

# Given materials
sodium = Material(780, 60, 1.7E-4, 1300)
flibe = Material(1940, 1, 2E-3, 2410)
fuel = Material(8530, 6, None, 500)

def fric_pressure_drop(g, mat, f = models.mcadams, lz = L, lc = DFLOW):
	"""Find the pressure drop due to friction over the channel
	
	Parameters:
	-----------
	g:          float, kg/s/m^2; mass flux
	mat:        Material; which material to use for the coolant
	f:          function; which model to use for friction factor
				[Default: models.mcadams]
	lz:         float, m; length of the channel in which coolant flows
				[Default: L]
	lc:         float; characteristic length for the flow
				[Default: DFLOW]
	
	Returns:
	--------
	float, Pa; pressure drop across
	"""
	re = mat.get_reynolds(g, lc)
	fff = f(re)
	return fff*lz/lc*g**2/(2*mat.rho)

print("\nPart 1: Friction pressure drop")

dpna = lambda g: fric_pressure_drop(g, sodium) - DPFRIC
gna = fsolve(dpna, 10)[0]   # Mass flux, sodium
dpms = lambda g: fric_pressure_drop(g, flibe, models.blasius) - DPFRIC
gms = fsolve(dpms, 10)[0]   # Mass flux, flibe
print("\tMass fluxes (kg/s/m^2)")
print("\t\tNa: {:.2f},  \t\tMS: {:.2f}".format(gna, gms))
rena = sodium.get_reynolds(gna)
rems = flibe.get_reynolds(gms)
print("\tReynolds numbers")
print("\t\tNa: {:.3e}, \t\tMS: {:.3e}".format(rena, rems))
mna = gna*AFLOW     # mass flow rate, sodium
mms = gms*AFLOW     # mass flow rate, flibe
print("\tMass flow rates (kg/s)")
print("\t\tNa: {:.3f},  \t\tMS: {:.3f}".format(mna, mms))

print("\nPart 2: Pump work")
wna = mna/sodium.rho*DPFRIC
wms = mms/flibe.rho*DPFRIC
print("\tNa: {:.2f} W, \tMS: {:.2f} W".format(wna, wms))

print("\nPart 3: Max Linear Power rate")
# Find the dimensions of equivalent annular cell
area = 1/2*sqrt(3)*DFLAT**2     # Area of a regular hexagon
d_eq = 2*sqrt(area/pi)          # diameter of a circle
print("\tEquivalent diameter: {:.2f} cm".format(d_eq*100))
dratio = (d_eq/DFLOW)
fcoeff = (2*d_eq**2/(d_eq**2 - DFLOW**2)*log(d_eq/DFLOW) - 1)
print("\tPrandtl number")
print("\t\tNa: {:.2e}, \tMS: {:.2f}".format(sodium.pr, flibe.pr))

def linear_power(delta_t, mat, nu, mdot):
	"""Find the linear heat generation rate to produce a desired
	temperature change, using the thermal resistors metaphor.

	Parameters:
	-----------
	delta_t:    float, K; temperature change in the fluid
	mat:        Material; which material to use for the coolant
	nu:         float; Nusselt number to use for heat transfer
	mdot:       float, kg/s; mass flow rate of the coolant

	Returns:
	--------
	q1:         float, kW/m; linear heat generation rate
	"""
	global fcoeff, fuel
	r_fuel = fcoeff/(4*pi*fuel.k)
	r_wall = 1/(pi*nu*mat.k)
	r_flow = L/(mdot*mat.c)
	q1 = delta_t / (r_fuel + r_wall + r_flow)
	return q1

dt = TMAX - TIN
nuna = models.lyon(rena, sodium.pr)
nums = models.dittus_boelter(rems, flibe.pr)
print("\tNusselt number")
print("\t\tNa: {:.2f}, \tMS: {:.2f}".format(nuna, nums))
q1na = linear_power(dt, sodium, nuna, mna)
q1ms = linear_power(dt, flibe, nums, mms)
print("\tMaximum q' (kW/m)")
print("\t\tNa: {:.2f}, \tMS: {:.2f}".format(q1na/1000, q1ms/1000))

print("""\nPart 4: Coolant selection

Flibe is the clear choice of coolant for the AHTR. It offers a {:.0%} advantage
in linear heat generation rate over sodium, while requiring a pump of only
{:.1%} of the capacity required by sodium.""".format(q1ms/q1na, wms/wna))
