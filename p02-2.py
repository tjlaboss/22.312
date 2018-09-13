# 22.312, PSet02
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 2-2: Relationships between assemblies of different pin arrays

"""
A utility wishes to replace the fuel in its existing PWR from 15 × 15 fuel
pin array assemblies to 17 × 17 fuel pin array assemblies. What is the ratio
of the core average linear power, q′, in the new core to the old core, assum-
ing that reactor power, length, and number of fuel assemblies are main-
tained constant? Repeat for core average heat flux, q′′.
"""

import reactors
from copy import copy

pwr17 = copy(reactors.pwr)
pwr17.name = "PWR-17"

pwr15 = reactors.ReactorType("PWR-15", ppitch = 1.26, npins = 15)
pwr15.npins = 15

pin_ratio = (pwr17.npins/pwr15.npins)

# Active fuel length doesn't change.
# There are more pins, each one generating less heat.
# Therefore, q'(17)/q'(15) will be inversely proportional to (17**2/15**2)
q1frac = round(1/pin_ratio**2, 3)
print(" q'(17)/q'(15)  =", q1frac)

# Heat flux is proportional to the total surface area: pi*D*L*Npins
# Total power, pi, and L are the same, but there are more pins:
#   Npins(17)/Npins(15) = pin_ratio**2
# and each one is smaller:
#   D(17)/D(15) = 1/pin_ratio
# Therefore, q''(17)/q''(15) will be inversely proportional to (17**2/15**2)*(15/17)
q2frac = round(1/pin_ratio, 3)
print("q''(17)/q''(15) =", q2frac)



