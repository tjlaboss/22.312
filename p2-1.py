# 22.312, PSet02
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 2-1: Relations among fuel element thermal parameters in various power reactors

from reactors import *

# Compute the core average values of the volumetric energy-generation rate in the
# fuel (q′′′) and outside surface heat flux (q′′) for the 6 reactor types


reactor_list = (bwr, pwr, candu, htgr, agr, lmfbr)
for rx in reactor_list:
	print('\n' + rx.name)
	
	q3 = rx.get_volumetric_heat()
	q3_print = str( round(q3*100**3 / 1000, 1) )
	print("q''':", q3_print, "MW/m^3")
	
	q2 = rx.get_heat_flux()
	q2_print = str( round(q2*100**2, 1) )
	print("q'': ", q2_print, "kW/m^2")