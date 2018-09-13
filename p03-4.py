# 22.312, PSet02
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 3-4: Decay heat energy

"""
Using one of the following equations:

	P/P0 = 0.066 * [ (T - T_s)^{-.2} - T^{-.2} ]

	P/P0 = 0.066 * [ t_s^{-.2} - (t_s + TAU_s)^{-.2} ]

Evaluate the energy generated in a 3411 MWt PWR after
it has operated for 1 year at 75% power.
"""

import scipy.integrate
import reactors
pwr = reactors.pwr
pwr.power = 0.75*3411

# Variables/constants
TAU_s = 365*24*3600   # time of shutdown
t_s = (3600, 3600*24, 3600*24*30.4) # time since shutdown
labels = ("1 hour: ", "1 day:  ", "1 month:")

# Equation 3.70d
P = lambda t: 0.066*(t**(-0.2) - (t + TAU_s)**(-0.2))   # 't' is time after shutdown

for i in range(len(t_s)):
	frac = scipy.integrate.quad(P, 0, t_s[i])[0]
	power = pwr.power * frac / 1E6
	print(labels[i], round(power, 3), "TJ")
	
print("\nThe values would be higher using equation 3.71 because it includes the")
print("decay heat of actinides U239 and Np239.")
