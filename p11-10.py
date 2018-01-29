# 22.312, PSet08
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 11-10: Flow Dynamics of Nanofluids

import models
from math import pi


# Given constants
D = 0.025           # m;      diameter of flow channel
AFLOW = pi/4*D**2   # m^2;    flow area of channel
VDOT = 4E-4         # m^3/s;  volume flow rate
BETA = 0.05         #         volume fraction of nanoparticles
RHO_W = 1000        # kg/m^3; density of water
MU = 0.001          # Pa-s;   viscosity of water
RHO_A = 4000        # kg/m^3; density of alumina
C = 780             # J/kg-K; specific heat capacity of alumina

print("\nPart 1) Mass flow rates")
m_n = BETA*RHO_A*VDOT
m_w = (1 - BETA)*RHO_W*VDOT
m_tot = m_n + m_w
print("\t -> water: {},\tal: {},\ttotal: {}  kg/s".format(m_w, m_n, m_tot))

print("\nPart 2) Pressure gradient")
rho_m = BETA*RHO_A + (1 - BETA)*RHO_W
g = m_tot/AFLOW
re = g*D/MU
fff = models.blasius(re)

print("""\
	rho_m: {:.0f} kg/m^3
	Aflow: {:.2e} m^2
	G:     {:.0f} kg/s/m^2
	Re:    {:.2e}
	f:     {:.3f}
""".format(rho_m, AFLOW, g, re, fff))

def dp_grav(z):
	return rho_m*9.81*z

def dp_fric(z):
	return fff*z/D*g**2/(2*rho_m)

print("\t -> DP = ({:.0f})*z + ({:.0f})*z".format(dp_grav(1), dp_fric(1)))
print("\t         ^fric       ^grav")

