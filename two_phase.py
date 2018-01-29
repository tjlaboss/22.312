# Two-phase
#
# Module with commmon equations for two-phase flow


def alpha(x, rhof, rhog, s=1):
	"""Find the void fraction given the quality
	
	Parameters:
	-----------
	x:          float; quality. must be on [0, 1]
	rhof:       float, kg/m^3; liquid density
	rhog:       float, kg/m^3; vapor density
	s:          float; slip ratio
				[Default: 1 -- HEM approximation]
	
	Returns:
	--------
	float; void fraction
	"""
	assert 0 <= x <= 1, "Quality is not on [0, 1]!"
	return 1/(1 + s*rhog/rhof*(1 - x)/x)


def quality(a, rhof, rhog):
	"""Find the quality given the void fraction
	
	Parameters:
	-----------
	a:          float; void fraction. must be on [0, 1]
	rhof:       float, kg/m^3; liquid density
	rhog:       float, kg/m^3; vapor density
	
	Returns:
	--------
	float; quality
	"""
	return 1/(1 + rhof/rhog*(1-a)/a)


def mixture_density(rhof, rhog, a):
	"""Density of a two-phase mixture
	
	Parameters:
	-----------
	rhof:       float, mass/volume; density of phase 1 (e.g., liquid)
	rhog:       float, mass/volume; density of phase 2 (e.g., vapor)
	a:          float; void fraction. Must be on [0, 1].
	
	Returns:
	--------
	float, mass/volume; density of the mixture in the same units
	as the densities given
	"""
	assert 0 <= a <= 1, "Void fraction is not on [0, 1]!"
	return a*rhog + (1 - a)*rhof
