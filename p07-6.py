# 22.312, PSet05
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 7-6: Containment sizing for a gas-cooled reactor with passive emergency cooling

from scipy.optimize import fsolve
import scipy.integrate as integrate

# Initial parameters
Q0 = 300        # MWt; reactor thermal power
Qcool = 0.02*Q0 # MWt; reactor cooling rate
P_MIN = 1.3E6   # Pa

# Thermodynamic parameters
# Primary system
V_p = 200   # m^3
T_p0 = 673  # K
P_p0 = 7E6  # Pa
# Conatinment
P_c0 = 1E5  # Pa (just below 1 atm)
T_c0 = 300  # K

# Helium properties. (Ideal gas)
HA = 0.004  # kg/mol
R = 8.31  # J/mol-K
c = 12.5  # J/mol-K

# Part 1: Find the containment volume so that the pressure does not fall below P_MIN
# immediately after a LBLOCA, assuming thermal equilibrium

"""
Unknowns:
	1. T2
	2. n_containment
	3. V_containment

Equations:
	1. Ideal gas law (intitial containment)
		P*Vc = nc*R*T0
	2. Ideal gas law (containment + primary)
		P*(Vp + Vc) = np*R*T2
	3. Conservation of Energy
		np*c*(T2 - T_p0) + nc*c*(T2 - T_c0) = 0
"""

# Ideal gas law: PV = nRT
# ...for primary
np = P_p0*V_p/(R*T_p0)
# ...for containment
V = lambda n_c: n_c*T_c0*R/P_c0  # volume of containment
# ...for entire volume (contaiment + primary)
T = lambda n_c: P_MIN*(V_p + V(n_c))/((np + n_c)*R)

# Conservation of energy
coe = lambda n_c: np*(T(n_c) - T_p0) + n_c*(T(n_c) - T_c0)
# Now everything is known!
nc = fsolve(coe, np)[0]
ntotal = nc + np
T1 = T(nc)
v_c = V(nc)

# Report
print("\nPart 1:")
print("\tT1  = {:.4} K".format(T1))
print("\tv_c = {:.4} m^3".format(v_c))


# Part 2: calculate at what time the pressure in the containment reaches its
# peak value after the LOCA as well as the peak temperature and pressure.
"""
Unknowns:
	1. time of maximum temperature
	2. maximum temperature
	3. maximum pressure

Equations:
	1. Reactor decay heat model
	2. Conservation of Energy
	3. Gay-Lussac's Law or Ideal Gas Law
"""



q = lambda t: 0.066*t**(-0.2)       # 't' is time after shutdown after infinite operation
dedt = lambda t: Q0*q(t) - Qcool    # Energy balance for the cooling rate (dE/dt)
tguess = 500                        # Seconds
tau = fsolve(dedt, 500)[0]          # Time at which the maximum temp. occurs

# Find the energy stored in the gas at tau
u_max = 1E6*integrate.quad(dedt, 0, tau)[0]     # Joules
# Then solve for temperature at an internal energy of u_max
# u_max = ntotal*c*(T2 - T1)
T2 = T1 + u_max/(ntotal*c)
# According to Gay-Lussac's law, the max pressure occurs at the max temperature like this:
P2 = P_MIN*(T2/T1)

# Report
print("\nPart 2:")
print("\ttau  = {} s".format(int(round(tau))))
print("\tTmax = {} K".format(int(round(T2))))
print("\tPmax = {:.3} MPa".format(round(P2*1E-6, 1)))


# Part 3
print("""
Part 3:
	To reduce the peak pressure in the containment, a nuclear engineer suggests venting the
containment gas to the atmosphere through a filter. What would be the advantages and
disadvantages of this approach?
	
Advantages:
	* The atmosphere would be an ultimate heat and pressure sink in a
		beyond-design-basis event
	* Containment could be made smaller and cheaper

Disadvantages:
	* The primary system already vents into the containment. Should fission gases
		be released in an accident, they would be vented to the atmosphere
	* If the containment is built smaller and the vents fail, pressure could build up
		beyond the design limit
""")
