# 22.312, PSet09
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 12-1

from iapws import IAPWS97 as Steam

# consts
PATM = 0.1013   # MPa
RE = 1E-5       # m
# Water
water = Steam(P=PATM, x=0)
vapor = Steam(P=PATM, x=1)
hfg = (vapor.h - water.h)*1000  # specific ethalpy, J/kg
dt_water = lambda r: 2*water.sigma*vapor.T*vapor.v/hfg/r
# Sodium.
# Proprties taken from: http://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
# Linearly interpolated at 1150 kelvins
tsat = 1151                              # K
sigma = 0.1203                           # N/m
hfg_na = (5079 + 5111 - 1273 - 1146)*500 # J/kg
rhov_na = (394 + 168.1)*0.0005           # kg/m^3
vg = 1/rhov_na                           # m^3/kg
dt_sodium = lambda r: 2*sigma*tsat*vg/hfg_na/r

print("\nPart 1) Relation between superheat and radius")
print("\tWater:  {:.2e}/re".format(dt_water(1)))
print("\tSodium: {:.2e}/re".format(dt_sodium(1)))

print("\nPart 2) Temperature and pressure difference with 10 um bubbles")
dtw = dt_water(RE)
dpw = hfg*dtw/(vapor.T*vapor.v)
print("\tWater:  dt(10 um) = {:.2f} K,\tdP = {:.1f} kPa".format(dtw, dpw/1000))
dtna = dt_sodium(RE)
dpna = hfg_na*dtna/(tsat*vg)
print("\tSodium: dt(10 um) = {:.1f} K,\tdP = {:.1f} kPa".format(dtna, dpna/1000))
