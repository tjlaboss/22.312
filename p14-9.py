# 22.312, PSet07
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 14.9

import models
from math import *
from scipy.optimize import fsolve

# Given constants
# B4C Control Rod
RHOC = 2500     # kg/m^3
CC = 950        # J/kg-K
KC = 35         # W/m-k
# Water
RHOW = 730      # kg/m^3
MU = 9E-5       # Pa-s
CW = 5500       # J/kg-K
KW = 0.56       # W/m-K
PRW = CW*MU/KW  # (dimensionless); Prandtl number
# Dimensions
DO = 0.026      # m; guide tube diameter
DC = 0.02       # m; control rod diameter
L = 3.7         # m; control rod length
F = 0.017       # (dimensionless); friction factor
K1 = 0.25       # loss coefficient of contraction
K2 = 1.0        # loss coefficient of expansion
# Conditions
TIN = 284       # degC
PIN = 15.5      # MPa
Q2AVG = 80000   # W/m^2
MDOTPRIME = 0.3 # kg/s

print("Part 1: Minimum flow rate through the guide tube")
"""Find the hydraulic diameter
        4*A_flow       4*pi/4*(DO**2 - DC**2)
D_H =  ----------  =  -------------------------
         P_wet             pi*(DO + DC)
"""
dh = (DO**2 - DC**2)/(DO + DC)
print("\tHydraulic diameter:  {:.2f} cm".format(dh*100))
aflow = pi/4*(DO**2 - DC**2)
print("\tFlow area:           {:.2f} cm^2".format(aflow*1E4))
com = lambda g: (K1 + K2 + 1.25*F*L/dh)*g**2/(2*RHOW) + (RHOW - RHOC)*9.81*L
gmin = fsolve(com, 1000)[0]
print("\tMass flux:           {:.1f} kg/s/m^2".format(gmin))
mdotmin = gmin*aflow
print("\tMass flow rate:      {:.2f} kg/s".format(mdotmin))

print("\nPart 2: Temperature in the control rod under reduced flow conditions")
print("\tPrandtl number:      {:.2f}".format(PRW))
g = MDOTPRIME/aflow
print("\tReduced mass flux:   {:.1f} kg/s/m^2".format(g))
re = g*dh/MU
print("\tReynolds number:     {:.0f}".format(re))
fff = models.mcadams(re)
print("\tFriction factor f:   {:.3f}".format(fff))
nu = models.dittus_boelter(re, PRW)
print("\tNusselt number:      {:.2f}".format(nu))
h = KW/dh*nu
print("\tHeat transfer coeff: {:.1f} kW/m^2-K".format(h/1000))
# Referring back to PSet05, Eq. 14.22b
zcrit = L/pi*atan(DC*h/(MDOTPRIME*CW)*L)
print("\tHeight at max. T:    {:.1f} cm".format(zcrit*100))
# And again from the previous PSet, Eq. 14.19
q2peak = pi/2*Q2AVG
q1peak = q2peak*pi*DC
tco = lambda z: TIN + q1peak*(cos(pi/L*z)/(pi*DC*h) +
                              (sin(pi/L*z) + 1)*L/(pi*MDOTPRIME*CW))
tmax = tco(zcrit)
print("\tMaximum crd T:       {:.1f} degC".format(tmax))

# Plotting
from pylab import *
zvals = linspace(-L/2, L/2)
yvals = [tco(zval) for zval in zvals]
plot(zvals, yvals)
xlim([-L/2, L/2])
ylim([TIN, tmax*1.02])
xlabel("z (m)")
ylabel("T ($^\circ$C)")
suptitle("Problem 14.9, Part 2")
title("Temperature profile across crd")
grid()
show()
