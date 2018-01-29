# 22.312, PSet05
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 7-9: Drain tank pressurization problem

from iapws import IAPWS97 as steam

# Constants
EPS = 0.02      # How close to converge
Y   = 0.2       # Underrelaxation factor
# Given parameters:
# Tank itself
P_TANK0 = 3     # MPa
P_BURST = 10    # MPa
P_SAT0  = 15.4  # MPa
V_TOTAL   = 12  # m^3
V_WATER0  = 2   # m^3
V_VAPOR0  = V_TOTAL - V_WATER0      # m^3
# Water in tank
water0 = steam(P = P_TANK0, x = 0)
vapor0 = steam(P = P_TANK0, x = 1)
M_l0 = V_WATER0/water0.v            # kg
M_v0 = V_VAPOR0/vapor0.v            # kg
M0 = M_v0 + M_l0                    # kg
X0 = M_v0/(M_l0 + M_v0)
U0 = M_l0*water0.u + M_v0*vapor0.u  # kJ
# Flow in
MDOT = 3                            # kg/s
liquid_in = steam(P = P_SAT0, x = 0)
# Air in tank (for part C only)
M_AIR = 11.93   # kg
R_AIR  = 268    # J/kg-K
CV_AIR = 719    # J/kg-K
V_AIR0 = V_VAPOR0
T_AIR0 = vapor0.T



# Part A: Setup
print("""
Part A:
	Control volume: Tank = 12 m^3
Unknowns:
	1. Time of burst 'tb'
	2. Mass of coolant at 'tb'
	3. Quality at 'tb'
	4. Final temperature (Part C only)

Equations:
	1. COM: integral (dm/dt) on t=[0,tb]  in control volume
		--> m(tb) = mdot*tb + m0
	2. COE: u_burst = u_initial + u(mdot.tb)
	3. Steam tables: u_burst = steam(P_burst, quality)
	4. Ideal gas law (Part C only)
""")

# Part B: Solve for the elapsed time to rupture


# Set up for iteration
def iterate(p_air2 = 0, tguess = 1000.0, xguess = 0.1, T2 = liquid_in.T):
	# Initialize some guesses before iterating
	last_t = 1200
	#v2 = V_WATER0/mguess
	V_gas2 = V_VAPOR0
	while abs(tguess - last_t)/last_t > EPS:
		last_t = tguess
		mflow = MDOT*last_t
		mguess = M0 + mflow
		
		if p_air2:
			# Partial pressure of air at T2
			p_air2 = T2*M_AIR*R_AIR/V_gas2 *1E-6    # MPa
			psat = P_BURST - p_air2                 # MPa
			de_air = M_AIR*CV_AIR*(T2-T_AIR0) *1E-3 # kJ
		else:
			psat = P_BURST
			de_air = 0
		
		liquid2 = steam(P = psat, x = 0)
		vapor2 = steam(P = psat, x = 1)
		T2 = liquid2.T  # K
		vmix2 = lambda x: liquid2.v*(1 - x) + vapor2.v*x
		umix2 = lambda x: liquid2.u*(1 - x) + vapor2.u*x
		u_fg = vapor2.u - liquid2.u  # kJ/kg
		
		# Control volume
		xguess = (V_TOTAL/mguess - liquid2.v)/(vapor2.v - liquid2.v)
		# Conservation of energy
		u2 = umix2(xguess)
		mflow = (U0 - de_air - M0*u2)/(u2 - liquid_in.u)
		# This calculation is unstable. Underrelax it!
		# (I'd like to thank Sterling for telling me to try this)
		tguess = Y*mflow/MDOT + (1-Y)*tguess
		
	
	return mguess, xguess, tguess
		
		
		

mb, xb, tb = iterate(False)
print("\nPart B: No air in system")
print("\tm {:.5} kg".format(mb))
print("\tx {:.5}".format(xb))
print("\tt {:.5} s".format(tb))

print("\nPart C:", M_AIR, "kg of air in system")
mc, xc, tc = iterate(True)
print("\tm {:.5} kg".format(mc))
print("\tx {:.5}".format(xc))
print("\tt {:.5} s".format(tc))

