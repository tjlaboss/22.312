# 22.312, PSet05
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 7-1: Containment pressure analysis

from iapws import IAPWS97 as steam

# Constants
V_C = 50970     # m^3
R_AIR = 286     # J/kg-K
M_AIR = 5.9E4   # kg
CV_AIR = 719    # J/kg-K
EPS = 0.02

# Starting conditions from Example 7.1
P0 = 0.523 # MPa
Pa0= 0.137 # MPa
T0 = 415.6 # Kelvins
X0 = 0.505
m_pw = 2.1E5 # kg
m_pg = 2.11E5# kg
mp = m_pw + m_pg
sat_liquid0 = steam(T=T0, x = 0)
sat_vapor0 = steam(T=T0, x = 1)
st0 = steam(T = T0, x = X0)
V_vapor0 = m_pg*sat_vapor0.v
e_p0 = m_pw*sat_liquid0.u + m_pg*sat_vapor0.u
#print("Partial pressure of st0 =", st0.P)

# Saturated liquid released from the secondary system
PS = 6.89 # MPa
VS = 89   # m^3
sts = steam(P = PS, x = 0)
print("Temperature of the secondary coolant at break: {:.4} K".format(sts.T))
ms = VS / sts.v
e_s0 = ms*sts.u
#print("Mass of secondary coolant released at {:.4} K: {:.4} kg".format(sts.T, ms))
m_w = ms + mp   # Total mass of coolant (primary + secondary)

# Control volume: All of the coolant in the system
V_TOTAL = V_C + VS + 354

# Set up the iteration
guess = steam(T = (sts.T + T0)/2, x = X0)
T2 = T0
count = 0
# While the current temperature is more than 1% away from the last guess:
while True:
	# Search for the right temperature
	# Set up saturated liquid and vapor phases for reference
	liquid = steam(T = guess.T, x = 0)
	vapor = steam(T = guess.T, x = 1)
	u_fg = vapor.u - liquid.u
	"""CONSERVATION OF ENERGY
	total_mass*(X*vapor.energy + (1-X)*liquid.energy)
	 - (primary0.energy + secondary0.energy)
	 - air.delta_energy = 0
	
	Solving for X:
		        primary_liquid.u0 + primary_vapor.u0 - air.delta_energy
		X = ---------------------------------------------------------------------
		                  total_coolant_mass * u_fg2
	"""
	de_air = M_AIR*CV_AIR*(guess.T - T0)/1000
	xguess = (e_p0 + e_s0 - de_air) / (m_w * u_fg) - liquid.u / u_fg
	# Check if the specific volume matches up
	vg = liquid.v + (V_vapor0 / m_w - liquid.v) / xguess
	
	T2 = guess.T
	if vg > (1 + EPS)*vapor.v:
		# Temperature is too high
		T2 *= 1 - EPS/2
	elif vg < (1 - EPS)*vapor.v:
		# Temperature is too low
		T2 *= 1 + EPS/2
	else:
		print("Converged after {} iterations.".format(count))
		break
	
	try:
		guess = steam(T = T2, x = xguess)
	except NotImplementedError:
		import sys
		print("\n\n**Error: X is too large; steam would be superheated.")
		break
		
	count += 1
	if count >= 1000:
		print("Did not converge after 1000 iterations; stopping")
		break


# Now account for the air pressure
liquid2 = steam(T = guess.T, x = 0)         # Saturated liquid at T2
V_liquid2 = m_w*(1 - guess.x)*liquid2.v     # Volume of said liquid
V_vapor2 = V_TOTAL - V_liquid2              # Volume of all vapor/gas
p_a2 = M_AIR * R_AIR * guess.T / V_vapor2 * 1E-6    # partial pressure of air
p2 = guess.P + p_a2                         #

print("\n7.1\tPressure({:.4} K): {:.4} MPa. Quality: {}%".format(guess.T, p2, round(100*guess.x)))
print("""\n7.2\tYes, it is possible for the containment pressure to drop, because the
\tsecondary coolant may act as a heat sink.
""")


