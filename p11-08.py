# 22.312, PSet08
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 11-08: Void, quality, and pressure drop problem

from iapws import IAPWS97 as Steam
import models
import two_phase
from math import pi

# Given constants
L = 3
D = 0.01
P = 7.4
MDOT = 0.3
MU = 8.7E-5
DBUBBLES = 0.0025
NBUBBLES = 4
VBUBBLES = NBUBBLES*pi/6 *DBUBBLES**3
AFLOW = pi/4*D**2
G = MDOT/AFLOW
RE = G*D/MU
# Steam properties
water = Steam(P=P, x=0)
vapor = Steam(P=P, x=1)

# Part 1
print("Flow area: {:.2e} m^2 \t Mass flux: {:.0f} kg/s/m^2".format(AFLOW, G))
print("Reynolds number:   {:.2e}".format(RE))
print("Bubble volume:     {:.2e} m^3".format(VBUBBLES))
vflow = 0.01*AFLOW
beta = VBUBBLES/vflow
alpha1 = beta
x = two_phase.quality(alpha1, water.rho, vapor.rho)
# Part 3
rho = two_phase.mixture_density(water.rho, vapor.rho, alpha1)
print("rho (mixture):     {:.1f} kg/m^3".format(rho))
fff = models.mcadams(RE)
print("friction factor:   {:.3f}".format(fff))
delta_p = rho*9.81*L + fff*L/D*G**2/(2*rho)

print()
print("1) Void fraction:  {:.3f}".format(alpha1))
print("   Quality:        {:.5f}".format(x))
print("3) Pressure drop:  {:.1f} kPa".format(delta_p*1E-3))