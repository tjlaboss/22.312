# 22.312, PSet08
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 5-1: Area-averaged parameters

import two_phase
from iapws import IAPWS97 as Steam

# Given constants
P = 7.2     # MPa
X = 0.15
S = 1.5
A = 0.012   # m^2
MDOT = 17.5 # kg/s
# Steam tables
water = Steam(P=P, x=0)
rhof = water.rho
vapor = Steam(P=P, x=1)
rhog = vapor.rho

alpha = two_phase.alpha(X, rhof, rhog, s=S)
beta = two_phase.alpha(X, rhof, rhog, s=1.0)
mdotv = X*MDOT
mdotl = (1-X)*MDOT
jayv = mdotv/(rhog*A)
gv = mdotv/A
gl = mdotl/A

print("""
alpha = {alpha:.4f}     mdotv = {mdotv:.2f} kg/s
beta  = {beta:1.4f}     mdotl = {mdotl:.2f} kg/s
jayv  = {jayv:.2f} m/s   Gv    = {gv:.2f} kg/s/m^2
                   Gl    = {gl:.2f} kg/s/m^2
""".format(**locals()))