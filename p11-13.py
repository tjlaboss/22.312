# 22.312, PSet08
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 11-13: SBLOCA in a BWR

from math import sqrt
from iapws import IAPWS97 as Steam

# Unit conversions
KJperKG_TO_BTUperLB = 0.43      # KJ/kg   -> BTU/lbm
MPA_TO_PSI = 145                # MPA     -> psi
LBperFT2_TO_KG_per_M2 = 4.882   # lb/ft^2 -> kg/m^2
# Given constants
A = 0.001   # m^2; area of the break
P0 = 0.1    # MPa; Containment pressure
P1 = 6.9    # MPa; primary coolant pressure
water1 = Steam(P=P1, x=0)   # primary coolant

def g(p):
	"""Mass flux as a function of pressure
	
	Parameters:
	-----------
	p:          float, MPa; critical pressure to find the mass flux at.
				If not available, use containment pressure..
	
	Returns:
	--------
	float, kg/s/m^2; mass flux
	"""
	return 0.61*sqrt(2*rhof*(P1 - p))*1E3


# Part 1: Non-equilibrium model with no nozzle
rhof = water1.rho
g1 = g(P0)
mdot1 = g1*A
print("1) mdot = {:.1f} kg/s".format(mdot1))

# Part 2: Non-equilibrium model with a short nozzle
fcr = 0.26  # based on a flimsy reading of Fig. 11.23
g2 = g(fcr*P1)
mdot2 = g2*A
print("2) mdot = {:.1f} kg/s".format(mdot2))

# Part 3: Equilibrium model with Moody's slip ratio
print("3) equilibrium model...")
p1eng = P1*MPA_TO_PSI
print("\tP1: {:.0f} psi".format(p1eng))
h0si = water1.h
h0eng = h0si*KJperKG_TO_BTUperLB
print("\th0: {:.1f} kJ/kg = {:.1f} BTU/lbm".format(h0si, h0eng))
g3eng = 7800    # sloppy chart reading, Fig 11.21
mdot3 = g3eng*LBperFT2_TO_KG_per_M2*A
print("\t -> mdot = {:.1f} kg/s".format(mdot3))
