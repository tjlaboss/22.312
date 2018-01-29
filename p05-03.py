# 22.312, PSet08
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 5-3: Estimating phase velocity differential

# Given constants
DTOT = 0.01
DBUB = DTOT/3
L = 3.66
G = 2000
P = 7.4
X = 0.0693
# Given steam properties
vg = 0.0239
vf = 0.001381

alpha = 4.0/9
print("alpha:   {:.3f}".format(alpha))
jayv = X*G*vg
print("jayv:    {:.2f} m/s".format(jayv))
vvv = jayv/alpha
print("vvv:     {:.2f} m/s".format(vvv))
jayl = (1 - X)*G*vf
print("jayl:    {:.2f} m/s".format(jayl))
vll = jayl/(1 - alpha)
print("vll:     {:.2f} m/s".format(vll))
print("delta v: {:.2f} m/s".format(vvv - vll))
print("delta j: {:.2f} m/s".format(jayv - jayl))
