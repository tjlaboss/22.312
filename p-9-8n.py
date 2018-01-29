# 22.312, PSet07
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 9-9N
# Pump seizure in a single-phase loop

import models
from math import pi
from scipy.optimize import fsolve

# Given constants
RHOL = 992.2        # kg/m^3
CP = 4200           # J/kg-K
K = 0.63            # W/m-K
MU = 6.5E-4         # Pa-s
BETA = 3.74E-4      # K^-1
# Temperatures and heat
THOT = 60           # degC
TCOLD = 20          # degC
DT = THOT - TCOLD   # K (or degC)
# Dimensions
H = 3.5             # m; height of thermal center
S = 4               # m; length of one side
LTOT = S*4          # m; total length of loop
D = 0.03            # m; inner diameter of loop


def pump(mdot):
	a = 1000
	b = 1.3
	dp = a*(1 - mdot/b)
	return dp

# Flow area inside loop
aflow = pi/4*D**2   # m^2
print("Aflow:                     {:.2f} cm^2".format(aflow*1E4))
# Boussinesq's Approximation
# Assume the given density is the lower one. It shouldn't really matter.
rhoc = RHOL
coef = 1 - BETA*DT
print("Fractional density change: {:.1%}".format(coef))
rhoh = coef*rhoc
print("Hot density:               {:.1f} kg/m^3".format(rhoh))
drho = rhoc - rhoh
print("delta_rho:                 {:.1f} kg/m^3".format(drho))
# Use Conservation of Momentum to solve for the mass flow rate
re = lambda m: 4*m/(pi*D*MU)
f = lambda m: models.blasius(re(m))
com_inf = lambda m: drho*9.81*H  - f(m)*LTOT/D * m**2/(2*RHOL*aflow**2)
com_0 = lambda m: com_inf(m) + pump(m)
mdot0 = fsolve(com_0, 1.0)[0]
print("\nBefore pump seizure:")
print("\tSteady-state mdot:     {:.2f} kg/s".format(mdot0))
# After the pump seizure, the pump power is 0
print("After pump seizure:")
mdotinf = fsolve(com_inf, 0.1)[0]
print("\tAsymptotic   mdot:     {:.2f} kg/s".format(mdotinf))