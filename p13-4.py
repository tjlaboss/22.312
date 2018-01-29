# 22.312, PSet10
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 13-4: Thermal parameters in a heated channel in two-phase flow

import math
from iapws import IAPWS97 as Steam

# Given constants
KELVIN = 273.15
L = 3
A = 1.5E-4
D = 2*math.sqrt(A/math.pi)
MDOT = 0.29
G = MDOT/A
P = 7.2
XOUT = 0.15
# Steam tables
water = Steam(P=P, x=0)
tsat = water.T
hin = water.h
vapor = Steam(P=P, x=1)
hfg = vapor.h - water.h
hout = (1 - XOUT)*water.h + XOUT*vapor.h

print("D:        {:.2f} cm".format(100*D))
print("G:        {:.1f} kg/s/m^2".format(G))
print("hin:      {:.1f} kJ/kg".format(hin))
print("hout:     {:.1f} kJ/kg".format(hout))
print("hfg:      {:.1f} kJ/kg".format(hfg))


print("\nPart 1: Fluid Temperature")
print("\t -> T = {:.1f} K = {:.1f} degC".format(tsat, tsat - KELVIN))

print("\nPart 2: Wall Temperature")
q2 = G*D*(hout - hin)/(4000*L)
print("\tq'':    {:.3f} MW/m^2".format(q2))
# Jens and Lottes Correletion
delta_t = 25*(q2/math.exp(4*P/6.2))**0.25
print("\tDeltaT: {:.1f} degC".format(delta_t))
twall = tsat + delta_t
print("\t -> Twall = {:.1f} degC".format(twall - KELVIN))

print("\nPart 3: Critical Power Ratio")
# Groeneveld correlation constants
k1 = math.sqrt(0.008/D)
k5 = 1
print("\tk1:     {:.3f} \t k5: {}".format(k1, k5))
# The rest of the stuff was from the lookup table...
qguess = 1.327
qlut = 1.919
qcorrected = k1*k5*qlut
qcr = (qcorrected + qguess)/2
print("\tq''_cr approximate:  {:.3f} MW/m^2".format(qcr))
cpr = qcr/q2
print("\t -> CPR ~ {:.2}".format(cpr))
