# 22.312, PSet07
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 10.6
# PWR Steam Generator

import models
from math import pi, log, exp
from scipy.optimize import fsolve

# Given
NT = 5700           # Number of tubes
LAVG = 16.0         # m; average length of steam generator tube
DO = 0.019          # m; Outer diameter of a tube
THICK = 0.0012      # m; Wall thickness of a tube
DI = DO - 2*THICK   # m; Inner diameter of a tube
AT = pi/4*DI**2     # m^2; flow area of a tube
KT = 26             # W/m-K; conductivity of the tube wall
PSEC = 5.6          # MPa; secondary loop pressure
TSEC = 272 + 273.15 # K; secondary loop saturation temperature
QT = 820E6          # W; heat transfer rate from primary -> secondary
MDOT = 5100         # kg/s; primary flow rate
GDOT = MDOT/(NT*AT) # kg/s/m^2; primary mass flux
# Water properties at 300 C and 15 MPA
RHO = 726           # kg/m^3; density
C = 5700            # J/kg-K; specific heat
MU = 92E-6          # Pa-s; viscosity
K = 0.56            # W/m-K; conductivity


# Solve
pr = C*MU/K
print("Prandtl number:  {:.3f}".format(pr))
# Heated perimeter
ph = pi*DI          # m; per tube
pht = NT*ph         # m; total over all tubes
# Heat transfer coefficient
re = GDOT*DI/MU
print("Reynolds number: {:.3e}".format(re))
nu = models.dittus_boelter(re, pr, cooled = True)
print("Nusselt number:  {:.1f}".format(nu))
htc = K/DI*nu
print("Heat transfer coefficient:  {:.2f} W/m^2-K".format(htc))
hw = 2*KT/(DI*log(DO/DI))
print("Heat transfer through tube: {:.2f} W/m^2-K".format(hw))
# Coefficient for the bulk temperature equation
m = MDOT/NT             # kg/s; mass flow rate per tube
alpha = (htc*ph)/(m*C*(1 + htc/hw))
# Finally, find the temperatures
print()
delta_t = QT/(MDOT*C)   # K; change in temperature, inlet to outlet
exponential = lambda tpout: TSEC + (tpout + delta_t - TSEC)*exp(-alpha*LAVG) - tpout
tout = fsolve(exponential, TSEC)[0] - 273.15
tin_ = tout + delta_t
print("Inlet temperature:  {:.1f} degC".format(tin_))
print("Outlet temperature: {:.1f} degC".format(tout))
