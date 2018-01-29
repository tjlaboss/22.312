# 22.312, PSet08
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 11-09: Calculation of some pipe's diameter

import models
import two_phase
from math import pi
from iapws import IAPWS97 as Steam
from scipy.optimize import fsolve

# Given constants
P = 6.7         # MPa; pressure
TFEED = 500     # K; feedwater temperature
QDOT = 856E6    # W; heater energy
H = 4.6         # m
XOUT = 0.1      # [dimensionless]; steam quality
# Steam tables
feed = Steam(P=P, T=TFEED)  # feedwater
hin = feed.h
print("hin:    {:.0f} kJ/kg".format(hin))
water = Steam(P=P, x=0)     # sat. liquid
hf = water.h
print("hf:     {:.0f} kJ/kg".format(hf))
vapor = Steam(P=P, x=1)     # proper steam
hg = vapor.h
print("hg:     {:.0f} kJ/kg".format(hg))
mix = Steam(P=P, x=XOUT)    # sat. mix

# Find the mass flow rate with COE
mst = QDOT/1000/(hg - hin)
mdot = mst/XOUT
# water mixture in the downcomer
hdown = hin*XOUT + hf*(1 - XOUT)
water_down = Steam(P=P, h=hdown)
rho_down = water_down.rho
print("hdown:  {:.0f} kJ/kg".format(hdown))
print("rho:    {:.1f} kg/m^3".format(rho_down))
dp_down = rho_down*H*9.81
print("DPdown: {:.1f} kPa".format(dp_down/1000))
# water+steam mixture in the riser
rhog = vapor.rho
rhof = water.rho
alpha = two_phase.alpha(XOUT, rhof, rhog)
print("alpha:  {:.4f}".format(alpha))
rhom = two_phase.mixture_density(rhof, rhog, alpha)
print("rhom:   {:.1f} kg/m^3".format(rhom))
dp_up = rhom*9.81*H
print("DPup:   {:.1f} kPa".format(dp_up/1000))
# friction loss in the riser
mu = water.mu
g = lambda d: 4*mdot/(pi*d**2)
re = lambda d: g(d)*d/mu
fff = lambda d: models.mcadams(re(d))
# Conservation of momentum
com = lambda d: dp_down - dp_up - fff(d)*H/d*g(d)**2/(2*rhom)
diameter = fsolve(com, 0.1)[0]

print("\n" + "-"*12)
print("mst:    {:.2f} kg/s".format(mst))
print("d:      {:.3f} m".format(diameter))
