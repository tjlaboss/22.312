# 22.312, PSet10
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 13-6
# Boiling crisis on the vessel outer surface during a severe accident

from math import pi
from iapws import IAPWS97 as Steam

# Given constants
KELVIN = 273.15
T0 = 80 + KELVIN
P0 = 0.1013          # MPa; atmospheric pressure
R = 5.2             # m;   vessel inner radius
THICK = 0.2         # m;   gap thickness
AHT = 2*pi*R**2     # m^2; area for heat transfer
Q2 = 350            # kW/m^2;  heat flux through the vessel
MDOT = 300          # kg/s;    mass flow rate
# Steam tables
sub_water = Steam(P=P0, T=T0)
hin = sub_water.h
sat_water = Steam(P=P0, x=0)
hf = sat_water.h
sat_vapor = Steam(P=P0, x=1)
hfg = sat_vapor.h - sat_water.h


print("AHT:    {:.2f} m^2".format(AHT))
qdot = AHT*Q2
print("qdot:   {:.1f} MW".format(qdot*1E-3))
hout = hin + qdot/MDOT
print("hout:   {:.1f} kJ/kg".format(hout))
x = (hout - hf)/hfg
print("x:      {:.3f}".format(x))
