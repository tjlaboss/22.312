# 22.312, PSet08
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 11-05: Impact of slip model on predicted void fraction

from iapws import IAPWS97 as Steam
import two_phase

# Constants
SHEM = 1.0
MDOT = 0.29
AFLOW = 1.5E-4
X = 0.15
P = 7.2
MPA_TO_PSI = 145.038
water = Steam(P=P, x=0)
vapor = Steam(P=P, x=1)


# Part 1: HEM Model
alpha1 = two_phase.alpha(X, water.rho, vapor.rho, SHEM)
print("1) Void fraction, HEM:     {:.3f}".format(alpha1))

# Part 2: Bankoff's Correction
psi = P*MPA_TO_PSI
k = 0.71 + 1E-4*psi
alpha2 = k*alpha1
print("2) Void fraction, Bankoff: {:.3f}".format(alpha2))

# Part 3: Dix's Correlation
# drift velocity for churn flow
vvj = 1.53*(water.sigma*9.81*(water.rho - vapor.rho)/water.rho**2)**0.25
jay = MDOT/AFLOW*((1 - X)/water.rho + X/vapor.rho)
b = (vapor.rho/water.rho)**0.1
c = alpha1*(1 + (1/alpha1 - 1)**b)
alpha3 = alpha1/(c + vvj/jay)
print("3) Drift flux model, Dix:  {:.3f}".format(alpha3))
