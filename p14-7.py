# 22.312, PSet10
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 14-7: Thermal hydralic analysis of a pressure tube reactor

import models
import two_phase
import math
from iapws import IAPWS97 as Steam
#from scipy.special import j0
from scipy.integrate import quad

# Math stuff
#jpi = 2.405
#RPEAK = 2.32
KELVIN = 273.15
# Given constants
KFUEL = 23          # W/m-K
KG = 50             # W/m-K
KC = 0.59           # W/m-K
CPC = 5000          # J/kg-K
RHOF = 776.3        # kg/m^3
RHOG = 35.94        # kg/m^3
MUC = 101E-6        # Pa-s
L = 6.0             # m
PITCH = 0.0275      # m
DFUEL = 0.0128      # m
RCOOL = 0.0074      # m
MDOT = 0.7          # kg/s
G = 1.4/(math.pi*RCOOL**2)  # kg/s/m^2
qdot = 260          # kW/fuel hole (hot channel)
# Steam properties
P = 6.89            # MPa
TIN = 245 + KELVIN  # K
hin = 1062.3        # kJ/kg
hf = 1261.6         # kJ/kg
hfg = 1511.9        # kJ/kg
sat_vapor = Steam(P=P, x=1)


q1max = qdot/L*math.pi/2
def q1(z):
	"""Linear heat generation rate in the hot channel
	
	Parameter:
	----------
	z:          float, m; axial distance. z=0 is midplane
	
	Returns:
	--------
	q':         float, kW/m; linear generation rate
	"""
	return q1max*math.cos(math.pi/L*z)
	
print("q'max:            {:.1f} kW/m".format(q1max))


# Precalculations
afuel = math.pi/4*DFUEL**2
acool = math.pi/2*RCOOL**2
agraf = math.sqrt(3)/4*PITCH**2 - afuel - acool
print("Area of fuel:     {:.2e} m^2".format(afuel))
print("Area of graphite: {:.2e} m^2".format(agraf))
print("Area of coolant:  {:.2e} m^2".format(acool))
d_eq = 2*math.sqrt((afuel + agraf)/math.pi)
dh = 2*RCOOL**2/d_eq
print("D_eq:             {:.1f} mm".format(d_eq*1000))
print("DH:               {:.2f} mm".format(dh*1000))
print("G:                {:.0f} kg/s/m^2".format(G))

print("\nPart 3: outlet enthalpy")
hout = hin + qdot/MDOT
print("\t -> hout = {:.1f} kJ/kg".format(hout))

print("\nPart 2: outlet temperature")
xout = (hout - hf)/hfg
print("\txout =         {:.3f}".format(xout))
if 0 <= xout <= 1:
	tout = sat_vapor.T
else:
	errstr = "Domain error: fluid is not saturated (x={:.4f})"
	raise ValueError(errstr.format(xout))
print("\t -> tout = {:.1f} degC".format(tout - KELVIN))

print("\nPart 4: outlet void fraction")
alpha = two_phase.alpha(xout, RHOF, RHOG, s=1)
print("\t -> alpha = {:.3f}".format(alpha))

print("\nPart 5: Non-boiling length")
hfrac = (hf - hin)/(hout - hin)
zonb = L/math.pi*math.asin(2*hfrac - 1)
lonb = L/2 + zonb
print("\t -> zonb = {:.2f} m = {:.2f} m (from inlet)".format(zonb, lonb))

print("\nPart 6: centerline temperature at onset of boiling")
re = G*dh/MUC
print("\tRe:            {:.2e}".format(re))
pr = CPC*MUC/KC
print("\tPr:            {:.2f}".format(pr))
nu = models.dittus_boelter(re, pr)
print("\tNu:            {:.1f}".format(nu))
htc = KC/dh*nu
print("\thtc:           {:.2f} kW/m^2-K".format(htc/1000))
source = 1000/(MDOT*CPC)*quad(q1, -L/2, zonb)[0]
tcl = q1(zonb)*(1/(4*math.pi*KFUEL) + math.log(d_eq/DFUEL)/(2*math.pi*KG) + \
                1/(math.pi*d_eq*htc))*1000 + TIN + source
print("\t -> TCL:       {:.1f} degC".format(tcl - KELVIN))

print("\nPart 7: Pressure drop across the channel")
h7 = lambda z: hf + 1/MDOT*quad(q1, zonb, z)[0]
xe7 = lambda z: (h7(z) - hf)/hfg
a7 = lambda z: two_phase.alpha(xe7(z), RHOF, RHOG)
rho7 = lambda z: two_phase.mixture_density(RHOF, RHOG, a7(z))
rhomout = rho7(L/2)
print("\trhom(L/2):    {:.1f} kg/m^3".format(rhomout))

print("\tDeltaP_grav:")
dpgrav_liq = RHOF*9.81*lonb/1000
dpgrav_mix = 9.81*quad(rho7, zonb, L/2)[0]/1000
dpgrav_tot = dpgrav_liq + dpgrav_mix
print("\t\t[-L/2, zonb] = {:.1f} kPa".format(dpgrav_liq))
print("\t\t[zonb, +L/2] = {:.1f} kPa".format(dpgrav_mix))
print("\t\tTotal: {:.1f} kPa".format(dpgrav_tot))

dpacc = (1/rhomout - 1/RHOF)*G**2/1000
print("\tDeltaP_acc:      {:.1f} kPa".format(dpacc))

print("\tDeltaP_fric:")
re7 = G*2*RCOOL/MUC
print("\t\tRe:           {:.2e}".format(re7))
fff = models.mcadams(re7)
print("\t\tfff:          {:.3f}".format(fff))
coeff = fff*G**2/(4000*RHOF*RCOOL)
dpfric_liq = coeff*lonb
print("\t\t[-L/2, zonb] = {:.1f} kPa".format(dpfric_liq))
phi_squared = lambda z: 1 + (RHOF/RHOG - 1)*xe7(z)
dpfric_mix = coeff*quad(phi_squared, zonb, L/2)[0]
print("\t\t[zonb, +L/2] = {:.1f} kPa".format(dpfric_mix))
dpfric_tot = dpfric_liq + dpfric_mix
print("\t\tTotal: {:.1f} kPa".format(dpfric_tot))
delta_p = dpfric_tot + dpacc + dpgrav_tot
print(" -> DeltaP = {:.0f} kPa".format(delta_p))
