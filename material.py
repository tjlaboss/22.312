# Material
#
# Classes for materials to be used in the simulation

import random


def number_density(density, A, frac = 1):
	"""Calculate the number density in atoms/cm^3 for each nuclie in the list.
	
	Inputs:
		density:    float; mass density of the material in g/cm^3
		A:          float; atomic mass (may be approximated as the number A)
		fracs:      float; atom fraction of the nuclide
					[Default: 1]
	
	Output:
		n:          float; number density in g/cm^3
	"""
	n = density * Avogadro / A * frac
	return n


class Nuclide(object):
	"""A simple nuclide with a name, mass, number density, and
	cross section dictionary
	
	Attributes:
		name:       str; human-readable name of the nuclide ("U235")
		mass:       float; atomic mass
		xs_dict:    dictionary of {"xs_type":xs}, where 'xs' is a float
		            representing the microscopic cross section (1/cm)
	"""
	def __init__(self, name, mass = 0.0, xs_dict = None):
		self.name = name
		self.mass = mass
		self.xs_dict = xs_dict
	
	def __str__(self):
		return self.name + " @ " + str(self.mass) + " g/mol"
		

class Material(object):
	"""A material with mass density, and nuclear composition.
	
	Parameters:
		name:       str; human-readable name of the material ("graphite")
		density:    float; mass density in g/cm^3
		nuclides:   dictionary of {Nuclide.name:Nuclide}
	
	Other attributes:
		wts:        dictionary of {Nuclide.name:weight fraction}
		ats:        dictionary of {Nuclide.name:atom fraction}
		total_xs:   float; total macroscopic cross section (1/cm)
		a_avg:		float; average of all the atomic masses weighted
					by their atom fractions
	"""
	def __init__(self, name, density = None, nuclides = None):
		self.name = name
		self.density = density
		self.nuclides = nuclides
		self.wts = {}
		self.ats = {}
		self.total_xs = 0.0
		self.a_avg = 0.0
	
	
	def __str__(self):
		rep = self.name
		rep += " @ " + str(self.density) + ":\n"
		rep += " containing the following nuclides:\n"
		for n in self.nuclides:
			rep += n + ', '
		return rep
	
	
	def convert_at_to_wt(self):
		"""Convert atomic fraction to weight fraction for this material's isotopes"""
		total_at = sum(self.ats.values())  # should be 1
		total_wt = 0.0
		
		if total_at >= 0:
			# already in weight fraction
			return
		else:
			for n in self.nuclides:
				total_wt += self.wts[n] * self.nuclides[n].mass
			for n in self.nuclides:
				if n not in self.wts:
					self.wts[n] = self.ats[n] * self.nuclides[n].mass / total_wt
	
	def convert_wt_to_at(self):
		"""Convert weight fraction to atomic fraction for this material's isotopes"""
		total_at = 0.0
		total_wt = sum(self.wts.values())
		
		if total_wt <= 0:
			# already in atomic fraction
			return
		else:
			for n in self.nuclides:
				total_at += self.wts[n] / self.nuclides[n].mass
			for n in self.nuclides:
				self.ats[n] = self.wts[n] / self.nuclides[n].mass / total_at
	
	
	def get_average_atomic_mass(self):
		"""Get self.a_avg: the average atomic mass
		If the value exists, return it; otherwise, calculate it.
		
		Output:
			self.a_avg: float; average of all the atomic masses
						weighted by their atom fractions
		"""
		if not self.a_avg:
			for n in self.nuclides:
				self.a_avg += self.nuclides[n].mass * self.ats[n]
		return self.a_avg
	
	
	def get_reaction_type(self):
		"""Choose the microscopic cross section for an interaction from the
		cross section dicionaries of self.nuclides.
		
		Outputs:
			xs_type:            str; key in 'cross_sections' for the reaction type
			xs:                 float; microscopic cross section (1/cm) for the reaction
		"""
		xi = random.random()
		total = self.get_total_xs()
		# New, less-efficient code that gets you the reaction type
		ratio = 0
		for n in self.nuclides:
			nuclide = self.nuclides[n]
			for xs_type in nuclide.xs_dict:
				xs = self.get_macro_xs(n, xs_type)
				ratio += self.ats[n]*(xs / total)
				if ratio > xi:
					return xs_type
		
		# If the for loop ends without returning a value, something went wrong!
		errstr = "The loop terminated without selecting a value;\n"
		errstr += "did you supply a list of floating point numbers?"
		raise TypeError(errstr)
		
	
	def get_macro_xs(self, nkey, reaction):
		"""Calculate the macroscopic cross section for a nuclide and reaction
		
		Inputs:
			nkey:       str; name of the nuclide (key in self.nuclides)
			rtype:      str; name of the reaction (key in Nuclide.xs_dict)
		
		Output:
			macro_xs:   float; macroscopic xs (1/cm)
		"""
		nuclide = self.nuclides[nkey]
		at = self.ats[nkey]
		# if reaction not in nuclide.xs_dict: return 0
		micro_xs = nuclide.xs_dict[reaction] * 1E-24
		A = self.get_average_atomic_mass()
		n = number_density(self.density, A, frac = at)
		macro_xs = micro_xs * n
		return macro_xs
		
	
	def get_total_xs(self):
		"""Get self.total_xs: the total macroscopic cross section.
		If the value exists, return it; otherwise, calculate it.
		
		Output:
			self.total_xs:  float; the total macroscopic xs (1/cm)
		"""
		if not self.total_xs:
			for n in self.nuclides:
				for r in self.nuclides[n].xs_dict:
					self.total_xs += self.get_macro_xs(n, r)
		return self.total_xs
	
