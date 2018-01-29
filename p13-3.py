# 22.312, PSet09
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 13-3: Heat transfer problems for a BWR channel

import models
from pylab import *
from scipy.optimize import fsolve
from iapws import IAPWS97 as Steam

PLOT = True
KELVIN = 273.15
# Given constants
L = 3.588               # m
D = 0.0112              # m
A = pi/4*D**2      # m^2
P0 = 7.14               # MPa
TIN = 278.3 + KELVIN    # K
G = 1625                # kg/s/m^2
Q1MAX = 47240           # W/m
# Steam tables
water_in = Steam(P=P0, T=TIN)
sat_water = Steam(P=P0, x=0)
sat_vapor = Steam(P=P0, x=1)
tsat = sat_water.T

hin = water_in.h
hf = sat_water.h
hg = sat_vapor.h
hfg = (hg - hf)
dhsat = (hf - hin)*1000

print("hin: {:.0f}\thf: {:.0f}\thg: {:.0f}\thfg: {:.0f}  kJ/kg".format(hin, hf, hg, hfg))
print("deltah, sat.:   {:.0f} J/kg".format(dhsat))
print("Tsat:           {:.1f} K".format(tsat))

print("\nPart 1: Axial location where equilibrium quality is 0")
xe0 = (hin - hf)/hfg
print("\tXe(-L/2):       {:.3f}".format(xe0))
qdot = 2/pi*Q1MAX*L
print("\tqdot:           {:.0f} kW".format(qdot/1000))
mdot = G*A
print("\tmdot:           {:.3f} kg/s".format(mdot))

def xe(z):
	return xe0 + qdot/(2*mdot*hfg*1000) * (sin(pi/L*z) + 1)


ze = fsolve(xe, -L/4)[0]
print("\t -> ze:        {:.2f} m = {:.2f} m from inlet".format(ze, ze + L/2))

print("\nPart 2: Axial location of ONB")

q2max = Q1MAX/(D*pi)
print("\tq''max:        {:.0f} kW/m^2".format(q2max/1000))
re = G*D/sat_water.mu
print("\tReynolds:      {:.2e}".format(re))
pr = sat_water.Prandt
print("\tPrandtl:       {:.2f}".format(pr))
nu = models.dittus_boelter(re, pr)
print("\tNusselt:       {:.0f}".format(nu))
htc = nu*sat_water.k/D
print("\th:             {:.1f} kW/m^2-K".format(htc/1000))

def tbulk(z):
	return TIN + L/(1000*pi*mdot*sat_water.cp) * \
	             Q1MAX*(sin(pi/L*z) + 1)

def q2(z):
	return q2max*cos(pi/L*z)

def twall(z):
	return tbulk(z) + q2(z)/htc

def delta_t(z):
	return twall(z) - tsat

def thom_onb(z):
	"""Find the heat flux as a function of z"""
	return 1E6*exp(2*P0/8.7)/22.7**2 * delta_t(z)**2

# Use Thom's correlation for heat flux to solve for the axial location
# where the onset of nuclear boiling occurs
fonb = lambda z: thom_onb(z) - q2(z)
zonb = fsolve(fonb, ze)[0]
print("\t -> zonb:      {:.2f} m = {:.2f} m from inlet".format(zonb, zonb + L/2))

if PLOT:
	figure()
	zvals = linspace(-L/2, L/2)
	xevals = array([xe(x) for x in zvals])
	plot(xevals, zvals, label="X_{eq}")
	plot([xe0, 0.5], [ze, ze], "gray")
	text(xevals.max()/1.5, -1.0, "$X_e = 0$")
	title("Equilibrium quality")
	ylabel("z (m)")
	xlim([xe0, 1.0])
	grid()

	fig2, ax1 = subplots()
	tbulkvals = array([tbulk(x) for x in zvals])
	ax1.plot(tbulkvals, zvals, "blue", label="$T_{bulk}$")
	twallvals = array([twall(x) for x in zvals])
	ax1.plot(twallvals, zvals, "orange", label="$T_{wall}$")
	ax1.plot([500, tbulkvals.max()], [zonb, zonb], "gray")
	text(twallvals.min() + 50, zonb + 0.15, "Onset of nucleate boiling")
	ax1.plot([], [], "green", label="$\Delta T$")
	ax1.set_xlim([500, twallvals.max()+10])
	legend()
	
	ax2 = ax1.twiny()
	dtvals = array([delta_t(x) for x in zvals])
	ax2.plot(dtvals, zvals, "green")
	ax2.set_xlim([0, 3*max(dtvals)])
	
	ylabel("z (m)")
	ax1.set_xlabel("T (K)")
	ax2.set_xlabel("$\Delta$T (K)", color="green")
	ax2.tick_params('x', colors='green')
	ax1.tick_params('y', colors="black")
	ax1.grid()
	show()
