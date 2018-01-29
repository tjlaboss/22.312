# Models
#
# Commonly-used thermal hydraulic correlations

def mcadams(re):
	"""McAdams correlation for friction factor

	Parameters:
	-----------
	:param re: Reynolds number

	Returns:
	--------
	:return: friction factor
	"""
	# assert 3E4 <= re <= 1E6, "Reynolds number {:.3e}".format(re) + \
	#	" is outside the domain of the McAdams correlation."
	return 0.184/re**0.2


def blasius(re):
	"""Blasius correlation for friction factor

	Parameters:
	-----------
	:param re: Reynolds number

	Returns:
	--------
	:return: friction factor
	"""
	# assert 4E3 <= re <= 1E5, "Re = {:.3e}".format(re) + \
	#	" is outside the domain of the Blasius correlation."
	return 0.316/re**0.25


def dittus_boelter(re, pr, cooled = False):
	"""The Dittus-Boelter correlation for the Nusselt number
	(dimentionless heat transfer coefficient)

	Parameters:
	-----------
	re:         float; Reynolds number
	pr:         float; Prandtl number (should be high for this correlation)
	cooled:     Boolean; whether the channel is cooled (True) or heated (False)
				[Default: True]

	Returns
	-------
	nusselt:    float; Nusselt number
	"""
	if cooled:
		n = 0.3
	else:
		n = 0.4
	nusselt = 0.023 * re**0.8 * pr**n
	return nusselt


def lyon(re, pr):
	"""The Lyon correlation for the Nusselt number

	Parameters:
	-----------
	re:         float; Reynolds number
	pr:         float; Prandtl number (should be high for this correlation)

	Returns
	-------
	nusselt:    float; Nusselt number
	"""
	nusselt = 7 + 0.025 * (re*pr)**0.8
	return nusselt


def film_boiling(water, vapor, tw, d):
	"""The film boiling correlation for the Nusselt number
	
	Parameters:
	-----------
	water:          iapws instance for saturated liquid (x=0)
	vapor:          iapws instance for saturated vapor (x=1)
	tw:             float, K; wall/surface temperature
	d:              float, m; hydraulic diameter
	
	Returns:
	--------
	nusselt:        float; Nusselt number
	"""
	dt = tw - water.T               # K; Temperature difference
	hfg = (vapor.h - water.h)*1000  # J/kg; enthalpy of vaporization
	hfg1 = hfg + 680*vapor.cp*dt    # J/kg
	rhofg = water.rho - vapor.rho   # kg/m^3
	nusselt = 0.62 * ( (9.81*vapor.rho*rhofg*hfg1*vapor.k**3) /
	                   (vapor.mu*d*dt) )**0.25
	return nusselt
