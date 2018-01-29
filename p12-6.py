# 22.312, PSet09
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 12-6: Analysis of decay heat removal during a severe accident

import models
import two_phase
from math import pi
from iapws import IAPWS97 as Steam

KELVIN = 273.15
SIGMAB = 5.67037E-8 # W/m^2/K^4
# Given constants
Q0 = 3400           # MWt
Q2 = 200E3          # W/m^2
TF = 2000 + KELVIN  # K
TSO = 112 + KELVIN  # K
TSAT = 100 + KELVIN # K
TIME = 3*3600       # seconds
D = 4.8             # m
R = D/2.0           # m
EP = 0.5            #
KVES = 30           # W/m-K
THICK = 0.22        # m
C0 = 1
# Steam tables
water0 = Steam(T=TSAT, x=0)
vapor0 = Steam(T=TSAT, x=1)

print("\nPart 1: Three hours after meltdown")
qd = 0.066*Q0*TIME**-0.2
print("\tqd:             {:.1f} MW".format(qd))
axs = pi*R**2
print("\tAxs:            {:.2f} m^2".format(axs))
nu = models.film_boiling(water0, vapor0, TF, D)
print("\tNusselt number: {:.1f}".format(nu))
hconv = vapor0.k/D*nu
print("\tConvective htc: {:.2f} W/m^2-K".format(hconv))
hrad = EP*SIGMAB*(TF**4 - TSAT**4)/(TF - TSAT)
print("\tRadiative htc:  {:.2f} W/m^2-K".format(hrad))
htc = hconv + hrad
qboil = htc*axs*(TF - TSAT)*1E-6
print("\tqboil:          {:.3f} MW".format(qboil))
hcond = KVES/THICK
print("\tConductive htc: {:.2f} W/m^2-K".format(hcond))
aht = 2*pi*R**2
print("\tA_HT:           {:.2f} m^2".format(aht))
qcond = hcond*aht*(TF - TSO)*1E-6
print("\tqcond:          {:.2f} MW".format(qcond))
dedt = qd - qboil - qcond
print("\tdE/dt:          {:.2f} MW".format(dedt))
print("\t -> system is", end = " ")
if abs(dedt) < 0.01*qd:
	print("holding steady")
elif dedt > 0:
	print("heating up")
else:
	print("cooling down")

print("\nPart 2: Void Fraction")
vvj = 1.53*(9.81*water0.sigma*(water0.rho - vapor0.rho)/water0.rho**2)**0.25
print("\tvvj:            {:.3f} m/s".format(vvj))
hfg = (vapor0.h - water0.h)*1000
g = Q2/hfg
print("\tG:              {:.3f} kg/s/m^2".format(g))
jayv = g/vapor0.rho
print("\tjayv:           {:.3f} m/s".format(jayv))
alpha = jayv/(C0 + vvj)
print("\t -> alpha = {:.3f}".format(alpha))

print("\nPart 3: HEM?")
print("\t -> beta = {:.3f}".format(two_phase.alpha(1, water0.rho, vapor0.rho)))
print("\t -> not a good approach")
