# xplot
#
# plot the quality as function of temperature

from iapws import IAPWS97 as steam
from pylab import *

V_TOT = 2500        # m^3
M_TOT = 2.12E6      # kg
VOL = V_TOT/M_TOT   # specific volume
T0 = 373            # K

def x(temp):
	vf = steam(T=temp, x=0).v
	vg = steam(T=temp, x=1).v
	return (VOL - vf)/(vg - vf)

temperatures = linspace(T0, T0+110)
qualities = array([x(ti) for ti in temperatures])
for i, q in enumerate(qualities):
	if q == max(qualities):
		break
txmax = temperatures[i]
plot(temperatures, qualities, label="X")
plot([temperatures.min(), temperatures.max()], [q, q], 'gray')
plot([txmax, txmax], [qualities.min(), qualities.max()], 'gray')
text(txmax, 1.05*q, "$X_{max} = $" + "{:.3e}".format(q) + " @ {:.4} K".format(txmax), ha="center")
vmin = steam(T=txmax,x=0).v*M_TOT*(1-q)
print("Minimum water volume: {:.5} m^3".format(vmin))
grid()
legend()
show()
