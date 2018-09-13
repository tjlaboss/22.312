# 22.312, PSet07
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 9.5
# The evil, evil problem that deprived me of sleep in undergraduate

import models
from math import pi
from scipy.optimize import fsolve

# Given constants
DP = 7.45E5 # Pa; pressure drop across assembly
ZO = 0.40   # m; height of the orifice
Z1 = 3.60   # m; height of assemblies in zone 1
Z2 = ZO + Z1# m; height of assemblies in zone 2
N1 = 100    # number of assemblies in zone 1
N2 =  80    # number of assemblies in zone 2
MU = 2E-4   # N-s/m^2; coolant viscosity
RHO = 800   # kg/m^3; coolant density
KC = 0.5    # loss coeficient of contraction
KE = 1.0    # loss coeficient of expansion
MDOT = 17.5E6 / 3600 # kg/s; total mass flow rate of coolant


# Conservation of energy w/ conservation of mass:
m2ass = MDOT/(2*N2) # kg/s; mass flow rate per assembly, zone 2
m1ass = N2/N1*m2ass # kg/s; mass flow rate per assembly, zone 1
print("mass per ass., zone 1: {:.2f} kg/s".format(m1ass))
print("mass per ass., zone 2: {:.2f} kg/s".format(m2ass))

# mass flux, zone 2
g2 = lambda dh: m2ass/(pi/4*dh**2)
# Reynolds number, zone 2
re2 = lambda dh: 4*m2ass/(pi*dh*MU)
# friction factor correlation, zone 2
f2 = lambda dh: models.mcadams(re2(dh))
# Conservation of Momentum, zone 2
com2 = lambda dh: RHO*9.81*Z2 + f2(dh)*Z2/dh*g2(dh)**2/(2*RHO) - DP
# Solve the CoM equation for hydraulic diameter of assemblies in zone 2
dh2 = fsolve(com2, 0.01)[0]
#print("Reynolds Number 2:    {:.3e}".format(re2(dh2)))
print("Hydraulic Diameter 2: {:.2f} cm".format(dh2*100))
mf2 = g2(dh2)
mf1 = N2/N1*mf2
print("G1: {:.2f}\tG2: {:.2f} kg/s/m^2".format(mf1, mf2))

# Zone 1
re1 = 4*m1ass/(pi*dh2*MU)
f1 = models.mcadams(re1)
dp1 = RHO*9.81*Z1 + f1*Z1/dh2*mf1**2/(2*RHO)
print("Pressure drop across Zone 1: {:.2e} Pa".format(dp1))
# Orifice
dpo = DP - dp1  # Pressure drop across just the orifice
go = lambda d: m1ass/(pi/4*d**2) / 4    # mass flux PER ORIFICE
reo = lambda d: 4*m1ass/(pi*d*MU)
fo = lambda d: models.mcadams(reo(d))
como = lambda d: RHO*9.81*ZO + (KC + KE + fo(d)*ZO/d)*go(d)**2/(2*RHO) - dpo
do = fsolve(como, 0.01)[0]

print("\nOrifice diameter: {:.2f} cm".format(do*100))
