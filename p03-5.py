# 22.312, PSet02
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 3-5: Decay heat from a PWR fuel rod


import scipy.integrate
from scipy.optimize import fsolve

P = lambda t: 0.066*t**(-0.2)    # 't' is time after shutdown after infinite operation

# 3411 MW core, 193 Assemblies of 264 fuel pins each; 97.4% of heat deposited in fuel
q0_rod = .974 * 3411E3 / (193 * 264) # kW
q_cool = 1  # W

# Energy balance: qrod - q_cool
f = lambda t: q0_rod * P(t) - q_cool

# The integral of f is the the total energy in the rod
# The derivative of the integral of f, set to 0, solves for the time with max heat in rod
# Therefore, set f itself to 0 and solve for t
guess = 1000 # seconds
tau = fsolve(f, guess)[0]

print("t:   {0} s".format(round(tau, 1)))


# Then integrate 'f' from [0, tau] to find the max value of q_rod
#tau = 1490.6
q_max = scipy.integrate.quad(f, 0, tau)[0]
print("Q_max: {0} kJ".format(round(q_max)))