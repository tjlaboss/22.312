# 22.312, PSet05
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 7-14:
# Behavior of a fully contained pressurized pool reactor under decay power conditions

import math
from iapws import IAPWS97 as steam
import pylab

# Reactor parameters
Q0 = 1600E3                             # kWt; steady-state core power
Q_DOT = 25E3                            # kWt; constant decay heat rate
TIME = 7*24*3600                        # seconds; time for cooling
# Geometry
R_VESSEL = 6.5                          # m
AXS_VESSEL = math.pi*R_VESSEL**2        # m^2; cross-sectional area of vessel
ZCORE = 2.0                             # m
ZCYL = 14.5                             # m
ZSTEAM0 = 2.17                          # m
Z_ABOVE = ZCYL - ZSTEAM0 - ZCORE        # m
V_STEAM0 = AXS_VESSEL*ZSTEAM0           # m^3; volume of the initial steam
V_BELOW = math.pi*2/3*R_VESSEL**3       # m^3; volume below core
V_ZCORE = AXS_VESSEL*ZCORE              # m^3; volume around core
V_ABOVE = AXS_VESSEL*Z_ABOVE            # m^3; volume above core
V_LIQUID0 = V_BELOW + V_ZCORE + V_ABOVE # m^3; volume if the initial liquid water
V_VESSEL = V_LIQUID0 + V_STEAM0         # m^3; entire pressure vessel volume

# Thermal hydraulic parameters
P0 = 0.10135                    # MPa
liquid0 = steam(P = P0, x = 0)
vapor0  = steam(P = P0, x = 1)
M_LIQUID0 = V_LIQUID0/liquid0.v
M_VAPOR0  = V_STEAM0 / vapor0.v
M0 = M_LIQUID0 + M_VAPOR0
X0 = M_VAPOR0/M0
U0 = M_LIQUID0 + liquid0.u + M_VAPOR0*vapor0.u
u0 = liquid0.u*(1-X0) + vapor0.u*X0
v0 = V_VESSEL/M0

"""
Unknowns:
	1. T2 or P2  at one week
	2. X2        at one week
	3. U, V      at one week
	
Equations:
	1. Conservation of Energy:  Energy_0 + Q_in = Energy_Final
	2. Steam tables:            U, V = steam(T2 or P2, X2)
	3. Conservation of Mass:    m_vapor0 + m_liquid0 = M0 = m_vapor2 + m_liquid2
"""

"""
Conservation of energy:

M0*[liquid2.u*(1-X) + vapor2.u*X] = integral(QDOT, [0, TIME]) + U0
"""

def x_iterate(t, xguess = X0, Tguess = 400, eps = 0.02):
	"""Iterate to find x(t) and T(t)
	
	Inputs:
		t:          time in s
		xguess:     where to start x at
		Tguess:     where to start T at (K)
		eps:        fractionally, how close to get
					[Default: 0.02]
	
	Outputs:
		x:          quality
		T:          Temperature (K)
	"""
	# Q_DOT is given as constant; the integral is trivial
	Q_decay = Q_DOT*t
	u_target = u0 + Q_decay/M0
	
	
	#
	Tmin = int(liquid0.T)
	Tmax = Tmin + 305
	for Tguess in range(Tmin, Tmax, 5):
		#print("\n\n\nTguess:", Tguess)
		xguess = 0.5
		last_x = 10*xguess
		u1 = 0.0
		c = 0
		converged = True
		overflooded = False
		while abs(xguess - last_x)/xguess > eps or abs(u_target - u1)/u_target > eps:
			last_x = xguess
			liquid1 = steam(T = Tguess, x = 0)
			vapor1 =  steam(T = Tguess, x = 1)
			u_fg = vapor1.u - liquid1.u
			v_fg = vapor1.v - liquid1.v
			v1   = liquid1.v + v_fg*xguess
			
			# Revise the x guess based on specific volume
			xguess = (v0 - liquid1.v)/v_fg
			u1 = liquid1.u + u_fg*xguess
			#print("u_target {:.4}\tu_actual {:.4}".format(u_target, u1))
			#print("v_target {:.4}\tv_actual {:.4}".format(v0, v1))
			
			#xguess = fsolve(v, last_x)[0]
			#xguess = (u_target - liquid1.u)/u_fg
			#xguess = fsolve()
			#print("xguess: {:.4} / last_x {:.4}".format(xguess, last_x))
			#print()
			
			c += 1
			if c > 10:
				converged = False
				break
		
		if converged or v1 > v0 or v1 == 0:
			if converged:
				print("Converged")
			if v1 > v0 or v1 == 0:
				overflooded = True
				print("Overflooded at t = {:.3} hours!".format(t/3600))
			break
		
		
	
	return v1, overflooded


def get_z(volume):
	"""
	Input:
		volume:
	Outputs:
	    z           height of the water
	    covered:    whether the core is covered
	"""
	if volume > V_ABOVE:
		z = volume/AXS_VESSEL
		return z, True
	elif volume > V_BELOW:
		z = volume/AXS_VESSEL
		return  z, False
	else:
		errstr = "Height calculations are only enabled for cylinders at this time."
		raise ValueError(errstr)
	

n = 21
z_array = pylab.zeros((n, 1))
c_array = pylab.ones((n, 1))
z_covered = pylab.zeros((n, 1))
z_uncovered = pylab.zeros((n, 1))
times = pylab.linspace(0, TIME, n)
t_array = pylab.linspace(0, TIME/3600, n)
ztop = R_VESSEL + ZCORE
zmax = R_VESSEL + ZCYL
for i in range(n):
	#t_array[i] = times[i]/3600
	volume, flooded = x_iterate(t = times[i])
	volume *= M0
	if flooded:
		z_array[i:] = zmax
		z_covered[i:] = zmax
		break
	else:
		z, c = get_z(volume)
		z_array[i] = z
		if c:
			z_covered[i] = z
		else:
			z_uncovered[i] = z
			c_array[i] = 0

# Set up plot
pylab.plot(t_array, [zmax]*n, color='k', linestyle = "dashed", label = "top of vessel")
pylab.plot(t_array, z_covered, color='b', marker="o", label = "water level")
pylab.plot(t_array, z_uncovered, color='r', marker="x")
pylab.plot(t_array, [ztop]*n, color='g', label = "top of core")
pylab.legend(loc = "center right")
pylab.xlabel("t (hours)")
pylab.ylabel("z (m)")
pylab.xlim([0, TIME/3600])
pylab.ylim([R_VESSEL, zmax+1])
pylab.title("Height of water in core", {"fontsize":16, "fontweight":"bold"})
pylab.show()