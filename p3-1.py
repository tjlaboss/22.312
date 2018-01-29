# 22.312, PSet02
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 3-1: Thermal design parameters for a cylindrical fuel pin

""" 3.1.1

Equation we're trying to solve:

	              <q'''>
	<phi> =   ---------------
	           E_f * SIGMA_F
           
Rewrite this equation:

	                power / V_fuel
	<phi> =   ---------------------------
	           E_F * sigma_235_f * N_235
"""


import math
import reactors
import material

# Global constants
N_ASSM = 193
M_FUEL = 558.5E3*N_ASSM     # kg
MEV_FISSION = 190           # MeV/fission
JOULE_FISSION = 3.04E-11    # J/fission
MICRO_XS = 35E-24           # cm^2
AVOGADRO =  6.022140857e+23 # things/mole

# Reactor parameters
pwr = reactors.pwr
pwr.power = 3411E6 * 0.974   # 97.4% actually desposited in fuel
pwr.npins = 264
pwr.nassm = N_ASSM
pwr.total_pins = pwr.npins * pwr.nassm
pwr.height = 365.8
# NOTE: Although the problem specifies cladding diameter of 9.5 mm,
# the answers in the textbook are given using a diameter of 9.62 mm!
pwr.dco = 0.962 #0.95

# Fuel
u235 = material.Nuclide("U235", 235.0439)
u238 = material.Nuclide("U238", 238.0508)
uranium = material.Material("Uranium (3.25% enriched)")
uranium.nuclides = {"U235":u235, "U238":u238}
uranium.wts = {"U235":.0325, "U238":1-0.0325}
uranium.convert_wt_to_at()

o16 = material.Nuclide("O16", 15.9994)
uo2 = material.Material("UO2 (3.25% enriched)")
uo2.nuclides = {"U235":u235, "U238":u238, "O16":o16}
uo2.ats["O16"] = 2.0/3
uo2.ats["U235"] = uranium.ats["U235"]*1.0/3
uo2.ats["U238"] = uranium.ats["U238"]*1.0/3
uo2.convert_at_to_wt()


# Let's find the number density of U235
"""            N_AV
	N = rho * ------
	            A

But since we know the mass and want to factor out the volume anyway,
replace `rho` with M_FUEL/V_fuel, and multiply through by V_fuel
"""
A_tot = len(uo2.nuclides) * uo2.get_average_atomic_mass()
N_V_fuel = M_FUEL * AVOGADRO / A_tot
N_V_u235 = N_V_fuel * uranium.ats["U235"]

macro_f = N_V_u235 * MICRO_XS

# And now we can plug in
flux = pwr.power / (JOULE_FISSION * macro_f)
print("Flux: {0:E} n/cm^2/s".format(flux))


"""3.1.2: Evaluate the average volumetric heat generation rate in the fuel.

Assume that the fuel density is 90% of theoretical.
This means we'll need to scale the value of q''' by 0.9

	q''' = E_f * macro_f * flux * 0.9
"""

pwr.q3 = pwr.power / (pwr.total_pins * math.pi/4 * pwr.df**2 * pwr.height)
q3 = 0.9 * pwr.q3     # W/cm^3 to MW/m^3 is unity; density is 90% of theoretical
print("q''': {0} MW/m^3".format(round(q3, 3)))

"""3.1.3:  Calculate the average linear power of the fuel

This will not change appreciably with the decrease in density;
use the theoretical value of q'''.

q' = q''' * pi/4 * df^2
"""

pwr.q1 = pwr.q3 * math.pi/4 * pwr.df**2
q1 = pwr.q1 * 100 / 1E3  # convert from W/cm to kW/m
# I'm expecting around 17.6, the answer I'm getting is a factor of 10 lower
print("q':   {0} kW/m".format(round(q1, 3)))

"""3.1.4: Calculate the average heat flux in the cladding outer radius

This is defined by the cladding diameter, and shouldn't change based
on the density of the fuel if the power is known and cladding is fixed.

	q'' = q' / (dco * pi)   --> calculated by the Reactor class
"""
q2 = pwr.get_heat_flux() * 1000 / 100
print("q'':  {0} kW/m^2".format(round(q2, 3)))
