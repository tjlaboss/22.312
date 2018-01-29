# 22.312, PSet08
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 11-12: Analysis of a liquid-metal-cooled reactor vessel

from math import sqrt

# Given constants
PATM = 0.1E6        # Pa; atmospheric pressure
R_VESSEL = 2.5      # m; inner radius of the vessel & height of the bottom
THICK = 0.04        # m; thickness of the steel vessel
H0 = 15             # m; height of the liquid lead in the cylindrical part
H1 = R_VESSEL + H0  # m; total height of the liquid in the vessel
ACK = 0.001         # m^2; area of the crack
T0 = 673.15         # K; initial temperature
P0 = 0.5E6          # Pa; temperature at the top of the vessel
# Lead
TBOIL = 1750        # degC; boiling temperature
RHO_L = 10500       # kg/m^3; density
MU = 1.6E-3         # Pa-s; viscosity
# Steel
SMAX = 138          # MPa; max. allowable stress intensity factor
RHO_S = 7500        # kg/m^3; density
# Nitrogen
R = 297             # J/kg-K; gas constant for nitrogen
CP = 1039           # J/kg-K
CV = CP - R         # J/kg-K
GAMMA = 1.4

print("Part 1: Mass flow rate through a crack")
print("(1a) crack in the bottom of the vessel")
p1 = P0 + RHO_L*9.81*H1
print("\tPressure at the bottom:  {:.1f} MPa".format(p1/1E6))
v1 = sqrt(2*(p1 - PATM  )/RHO_L)
print("\tVelocity out the bottom: {:.1f} m/s".format(v1))
mdot1 = RHO_L*v1*ACK
print("\t -> mdot out the bottom: {:.1f} kg/s (lead)".format(mdot1))
print("(1b) crack at the top of the vessel")
fcr = (2/(GAMMA + 1))**(GAMMA/(GAMMA - 1))
print("\t{:.2f}x{:.1f} MPa > {:.1} MPa".format(fcr, P0/1E6, PATM/1E6))
rho0 = P0/(R*T0)
print("\trho0(nitrogen):          {:.2f} kg/m^3".format(rho0))
g2 = rho0*sqrt(2*CP*T0*( fcr**(2/GAMMA) - fcr**((GAMMA+1)/GAMMA) ) )
print("\tmass flux of nitrogen:   {:.0f} kg/s/m^2".format(g2))
mdot2 = g2*ACK
print("\t -> mdot out the top:    {:.4f} kg/s (nitrogen)".format(mdot2))
