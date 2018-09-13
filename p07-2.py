# 22.312, PSet05
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 7-2: Ice condenser containment analysis
# Calculate the mass of ice needed to stay below maximum containment pressure

from iapws import IAPWS97 as steam
from scipy.optimize import fsolve

# Fixed parameters
P_MAX   = 0.4       # MPa
V_C     = 5.05E4    # m^3 (containment volume)
C0      = 273.15    # 0 degrees Celcius
EPS     = 0.005     # Get within 0.5% of the answer
MAX_ITER= 100       # Prevent infinite loops

# Initial conditions
P0   = 0.1013       # MPa (containment pressure)
T_c0 = 300          # K   (containment temperature)
P_p0 = 15.5         # MPa (primary coolant pressure)
V_p0 = 500          # m^3 (primary coolant volume)
st0 = steam(P = P_p0, x = 0) # primary coolant, sat liquid
mp = V_p0/st0.v

# Ice
T_ice0 = 263        # K
lf_ice = 333        # kJ/kg
c_ice  = 4.230      # kJ/kg-K
c_water= 4.18       # kJ/kg-K
liquid_ice = steam(T = T_c0, x = 0)

# Air
R_AIR  = 268    # J/kg-K
CV_AIR = 719    # J/kg-K
# PV = mRT --> m = PV/RT
M_AIR = (P0*1E6 * V_C) / (R_AIR * T_c0)     # kg
MRV_AIR = M_AIR*R_AIR / (V_C - V_p0) * 1E-6 # MPa


"""Assumptions:
	* Dry containment (negligible relative humidity)
	* The entire primary loop is evacuated into containment
	* The entire block of ice melts
	* The final state will saturated
"""

"""
Unknowns:

1. mice -- what we're looking for
2. u2   -- required for energy balance
3. x    -- required for u2
4. PSat -- required for u2
5. T2   -- required for energy balance

Equations:
 
1. Conservation of energy:
	Q_loss(primary) = m*c*dT(air) + m_ice*[c_ice*(0 - Tice) + lf_melt + c_water(T0 - 0) + {u(T2) - u(T0)}]

2-4. Steam tables:
	u2 = f(steam(x, PSat))
	T2 = f(steam(x, PSat))
	x  = f(u2) --> codependent with u2, T2; requires iteration

5. Ideal gas law:
	psat = P_MAX - (m*R*T2/V)_air
	 --> codependent with T2 and mice; requires iteration
"""


# Initial guesses
last_x = 1
x = 0.75
p_air2 = P0
last_T = 500
c = 0
# Search for the right mass, quality, temperature, and pressure
while abs(x - last_x) / last_x > EPS:
	# Partial pressure of air at T2
	p_air2 = last_T * MRV_AIR
	psat = P_MAX - p_air2
	# Saturated liquid and vapor states at max pressure
	liquid2 = steam(P = psat, x = 0)
	u_l = liquid2.u*(1-x)
	vapor2 = steam(P = psat, x = 1)
	u_v = vapor2.u * x
	T2 = liquid2.T
	
	# Energy calculations (function of T2, X)
	Q_primary = mp*(u_l + u_v - st0.u)
	Q_air = M_AIR * CV_AIR * (T2 - T_c0) * 1E-3
	
	q_ice0 = c_ice * (C0 - T_ice0)
	q_melt = lf_ice
	q_ice1 = c_water * (T_c0 - C0)
	q_ice2 = (liquid2.u*(1 - x) + vapor2.u*x) - liquid_ice.u
	
	Q_all = lambda m: m*(q_ice0 + q_ice1 + q_ice2 + q_melt) + Q_air + Q_primary
	mice = fsolve(Q_all, 1E5)[0]
	
	# Guess a new quality
	last_x = x
	# Use the mass of the ice to find the specific volume
	v_target = (V_C + V_p0) / (mp + mice)
	# Then, use specific volume to intelligently estimate a new x
	v_actual = liquid2.v * (1 - x) + vapor2.v * x
	#print("vg_target: {:.3}\tvg_actual: {:.3}".format(v_target, v_actual))
	x = (liquid2.v - v_target) / (liquid2.v - vapor2.v)
	
	
	last_T = T2
	if c > MAX_ITER:
		print(MAX_ITER, "iterations reached; canceling.")
		break
	else:
		c += 1

print("\nConverged after", c, "iterations.\n")

print("X    = {:.3} %".format(100*x))
print("u_w2 = {:.4} kJ/kg".format(u_l + u_v))
print("T2   = {:.4} K".format(T2))
print("Psat = {:.3} MPa".format(psat))
print("mice = {:.2E} kg".format(mice))
