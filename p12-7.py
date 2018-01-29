# 22.312, PSet11
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 12-7: Void fraction and pressure drop in an isolated condenser

from pylab import *
from scipy.optimize import fsolve
from iapws import IAPWS97 as Steam
import models
import two_phase

KELVIN = 273.15
# Given constants
N = 200             # tubes
L = 12              # m
D = 0.03            # m
AXS = pi/4*D**2     # m^2
AHT = pi*D*L        # m^2
MDOT = 50           # kg/s
G = MDOT/(N*AXS)    # kg/s/m^2
TSAT = 280 + KELVIN # K
VF = 0.0013         # m^3/kg
RHOF = 1/VF         # kg/m^3
VG = 0.03           # m^3/kg
RHOG = 1/VG         # kg/m^3
HF = 1237           # kJ/kg
HG = 2780           # kJ/kg
HFG = HG - HF       # kJ/kg
MUF = 9.8E-5        # Pa-s
MUG = 1.9E-5        # Pa-s
KF = 0.574          # W/m-K
KG = 0.061          # W/m-K
# Other steam table data
sat_water = Steam(T=TSAT, x=0)
cpf = sat_water.cp
sat_vapor = Steam(T=TSAT, x=1)
psat = sat_vapor.P


print("\nPart 1: Isolation condenser heat removal rate")
qdot = MDOT*HFG
print("\t -> Qdot = {:.2f} MW".format(qdot/1000))

print("\nPart 2: Inner surface temperature on the tube")
q2 = qdot/AHT
print("\tq'':     {:.1f} MW/m^2".format(q2/1000))

def chato(dt):
	"""Chato's correlation for the heat transfer coefficient
	
	Parameter:
	----------
	dt:         float, K; temperature differential between Twall and Tsat
	
	Returns:
	--------
	hbar:       float, W/m^2-K; avg. heat transfer coefficient for condensation
	"""
	hfg_star = HFG + 0.68*cpf*dt
	numer = 9.81*(RHOF - RHOG)*RHOF*KF**3*hfg_star
	denom = D*MUF*dt
	return 0.555*(numer/denom)**0.25

f = lambda dt: chato(dt) - q2/dt
delta_t = fsolve(f, 60)[0]
print("\tDeltaT:  {:.1f} degC".format(delta_t))
twall = TSAT - delta_t
print("\t -> Twall = {:.1f} K = {:.1f} degC".format(twall, twall - KELVIN))

# Part 3
def x(z):
	"""Flow quality as a function of axial position.
	Uses a purely linear relationship.
	
	Parameter:
	----------
	z:          float, m; axial position from the inlet
	
	Returns:
	--------
	x:          float; quality of the steam
	"""
	return z/L

zvals = linspace(0, L)
xvals = array([x(z) for z in zvals])
avals = array([two_phase.alpha(x, RHOF, RHOG) for x in xvals])
plot(zvals, xvals, "b-", label = "Quality $X$")
plot(zvals, avals, "r-", label = "Void fraction $\\alpha$")
xlim([0, L])
ylim([0, 1.05])
xlabel("z (m)")
title("Problem 3: Axial profile")
legend(loc='lower right')
grid()

print("\nPart 4: Pressure drops across the tubes")
print("\tG:            {:.0f} kg/s/m^2".format(G))
re = G*D/MUF
print("\tRe:           {:.2e}".format(re))
fff = models.mcadams(re)
print("\tftp:          {:.3f}".format(fff))
rhobar = 2/(VG + VF)
print("\trhobar:       {:.0f} kg/m^3".format(rhobar))
dp_acc = G**2*(VG - VF)/1000
print("\tDeltaP_acc:  {:+.2f} kPa".format(dp_acc))
dp_fric = -fff*L/D*G**2/(2000*rhobar)
print("\tDeltaP_fric: {:+.2f} kPa".format(dp_fric))
dp_tot = dp_acc + dp_fric
print("\t -> DeltaP_tot = {:+.2f} kPa".format(dp_tot))
show()
