# 22.312, PSet09
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 12-4: Shell and tube horizontal evaporator

from scipy import pi, sqrt
from scipy.optimize import fsolve
import models

# Given constants
PATM = 0.1013           # MPa; atmospheric pressure
KELVIN = 273            # K;   0 degC to Kelv
PIN = .69               # MPa; inlet pressureins
TIN = 130 + KELVIN      # K;   inlet temperature
TOUT = 120 + KELVIN     # K;   outlet temperature
TAVG = (TIN + TOUT)/2   # K;   average temperature
DT = 10                 # K;   temperature difference
PSH = PATM              # MPa; pressure on the shell side
TSH = 100 + KELVIN      # K;   temperature on the shell side
NTUBES = 30             # tubes
CSF = 0.013             # <dimensionless>
TFLUX = 125 + KELVIN    # K;   temperature used for heat flux calcs
D = 0.01                # m;   diameter of one tube
A = pi/4*D**2           # m^2; flow area of one tube
ATOT = NTUBES*A         # m^2; total flow area
VW = 3                  # m/s; velocity of water
# Steam tables--Liquid
RHOL = 960              # kg/s
CPL = 4200              # J/kg-K
MUL = 3E-4              # Pa-s
KL = 0.68               # W/m-K
PRL = 1.75              # (Prandtl number)
SIGMA = 0.06            # N/m
# Steam tables--Vapor
RHOV = 0.60             # kg/s
CPV = 2.0               # J/kg-K
MUV = 1.3E-5            # Pa-s
KV = 0.025              # W/m-K
PRV = 0.97              # (Prandtl number)
HFG = 2.28E6            # J/kg
# Other


# Easy stuff
G = RHOL*VW
print("Area:        {:.2e} m^2;\t Mass flux G: {:.0f} kg/s/m^2".format(A, G))
MTOT = G*ATOT
MDOT = G*A
print("mdot(tube):  {:.2f} kg/s;\t\t mdot(tot.):  {:.2f} kg/s\n".format(MDOT, MTOT))
rel = G*D/MUL
print("Reynolds number:   {:.2e}".format(rel))
fff = models.blasius(rel)
print("friction factor f: {:.3f}".format(fff))
qdot = MTOT*CPL*DT
print("Qdot:              {:.3f} MW".format(qdot*1E-6))
nu = models.dittus_boelter(rel, PRL, cooled = True)
print("Nusselt number:    {:.2f}".format(nu))
h = KL/D*nu
print("htc:               {:.2f} kW/m^2-K".format(h/1000))
#q2 = h*25
#print("Avg. heat flux:    {:.1f} kW/m^2".format(q2/1000))

# Rohsenow Correlation
rhofg = RHOL - RHOV
rohsenow = lambda q2: MUL*HFG*sqrt(9.81/SIGMA*rhofg) * \
                      (CPL*(TAVG - TSH - q2/h)/(CSF*HFG*PRL)) - q2
heat_flux = fsolve(rohsenow, 4E4)[0]
print("q'':               {:.2e} W/m^2".format(heat_flux))
aht = qdot/heat_flux
print("A_HT:              {:.2f} m^2".format(aht))
print()
l = aht/(NTUBES*D*pi)
print("Part 1: length of the tubes  = {:.2f} m".format(l))
edot = qdot/HFG
print("Part 2: evaporation rate     = {:.3f} kg/s".format(edot))
print("Part 3: total mass flow rate = {:.2f} kg/s".format(MDOT))
dp_tu = fff*l/D*G**2/(2*RHOL)
print("Part 4: DeltaP in the tubes  = {:.2f} kPa".format(dp_tu/1000))
