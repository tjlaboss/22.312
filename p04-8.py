# 22.312, PSet02
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 4-8: Internal conservation equations for an extensive property

from scipy import *
import scipy.integrate as integrate

# Given
V_MAX = 2.0 # 2.0 m/s
R = 0.05    # m
velocity = lambda r: V_MAX * (1 - (r/R)**2)

# Part 1: Coolant flow rate in the tube
#
# Equation for volumetric flow rate: Area * velocity
# Integrate that equation over 'r' from 0 to R
integrand1 = lambda r: velocity(r) * (2*pi*r)
vdot = integrate.quad(integrand1, 0, R)[0]
print("Coolant flow rate: {0:.3E} m^3/s".format(vdot))


# Part 2: Average coolant velocity
#
# Find the velocity by integrating over the velocity profile,
# and dividing by the integral of the radii
integrand2 = lambda r: (2*pi*r)
v_avg = integrate.quad(integrand1, 0, R)[0] / integrate.quad(integrand2, 0, R)[0]
print("Average velocity:  {0:.3} m/s".format(v_avg))

