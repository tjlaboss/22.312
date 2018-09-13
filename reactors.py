# Reactors
#
# Module containing parameters for typical BWR, PWR, CANDU, HTGR, AGR, AND LMFBR specimens

from math import *

class ReactorType(object):
	"""Container for reactor data
	
	Required Parameters:
		name:       str; name of the reactor
		ppitch:     float; pin pitch (cm)
		npins:      int; number of pins across an assembly
		
	Optional Parameters:
		df:         float; diamater of the fuel cylinder (cm)
		dci:        float; diameter of the inside of the glad (cm)
		dco:        float; diameter of the outside of the glad (cm)
		height:     float; active length of fuel (cm)
		apitch:     float; assembly pitch (cm)
		nassm:      int; number of assemblies
		shape:      str; {"square", "hexagon", or "circle"}
	"""
	def __init__(self, name, ppitch, npins,
	             df = None, dci = None, dco = None, height = None,
	             apitch = None, nassm = None, shape = "square"):
		self.name = name
		self.ppitch = ppitch
		self.npins = npins
		
		self._rfo = None
		self._dfo = None
		self._rci = None
		self._dci = None
		self._rco = None
		self._dco = None
		
		self.df = df
		self.dci = dci
		self.dco = dco
		
		
		self.height = height
		
		self.nassm = nassm
		self.apitch = apitch
		self.shape = shape
		
		# Heat generation parameters that can be changed by functions
		self.power = None   # MW
		self.Q3 = None      # Power density core average
		self.q1 = None      # Linear heat rate (avg)
		self.q2 = None
		self.q3 = None
		self.q1max = None   # Linear heat rate (max)
		
	def __str__(self):
		return self.name
	
	@property
	def dco(self):
		return self._dco
	
	@property
	def rco(self):
		return self._rco
	
	@property
	def dci(self):
		return self._dci
	
	@property
	def rci(self):
		return self._rci
	
	@property
	def dfo(self):
		return self._dfo
	
	@property
	def rfo(self):
		return self._rfo
	
	@dco.setter
	def dco(self, dco):
		self._dco = dco
		if dco is not None:
			self._rco = dco/2.0
		else:
			self._rco = None
	
	@dci.setter
	def dci(self, dci):
		self._dci = dci
		if dci is not None:
			self._rci = dci/2.0
		else:
			self._rci = None
	
	@dfo.setter
	def dfo(self, dfo):
		self._dfo = dfo
		if dfo is not None:
			self._rfo = dfo/2.0
		else:
			self._rfo = None
	
	@rco.setter
	def rco(self, rco):
		self._rco = rco
		if rco is not None:
			self._dco = rco*2.0
		else:
			self._dco = None
	
	@rci.setter
	def rci(self, rci):
		self._rci = rci
		if rci is not None:
			self._dci = rci*2.0
		else:
			self._dci = None
	
	@rfo.setter
	def rfo(self, rfo):
		self._rfo = rfo
		if rfo is not None:
			self._dfo = rfo*2.0
		else:
			self._dfo = None


	def get_volumetric_heat(self):
		if not self.q3:
			if self.q1 and self.dco:
				self.q3 = 4/pi * self.q1 / self.df**2
			elif self.q2:
				if self.dco and self.df and self.q2:
					if self.shape in ("square", "triangle"):
						self.q3 = self.q2 * self.dco / self.df**2
					else:
						raise TypeError("Cannot calculate this shape yet")
			else:
				raise TypeError("Not enough information to calculate q'''.")
		return self.q3
		
	
	def get_heat_flux(self):
		if not self.q2:
			assert self.q1, "Linear heat rate (q1) needs to be set"
			assert self.height, "Active fuel height (height) needs to be set"
			assert self.dco, "Clad outer diameter (dco) needs to be set"
			
			self.q2 = self.q1 / (self.dco * pi)
		return self.q2
	


# Common reactor types

# GE Boiling Water Reactor
bwr = ReactorType("BWR/5", ppitch = 1.437, npins = 9, df = 0.96, dco = 1.12,
                  apitch = 15.2, height = 358.8)
#bwr.power = 3323     # MW
bwr.Q3 = 54.1       # kW / L
bwr.q1 = 0.176      # kW / cm
bwr.q1max = 0.4724  # kW / cm

# Westinghouse Pressurized Water Reactor
pwr = ReactorType("PWR", ppitch = 1.26, npins = 17, df = 0.8192, dco = 0.95,
                  apitch = 21.5, height = 365.8)
#pwr.power = 3411    # kW
pwr.Q3 = 105.0      # kW / L
pwr.q1 = 0.1786     # kW / cm
pwr.q1max = 0.4462  # kW / cm

# Pressurized Heavy Water Reactor
candu = ReactorType("CANDU", ppitch = 1.46, npins = 37, df = 1.22, dco = 1.31,
                    shape = "circle", apitch = 28.6, height = 594)
#candu.power =
candu.Q3 = 12.0
candu.q1 = .257
candu.q1max = .441

# High Temperature Gas Reactor
"""Note:
Approximating the fuel as a cylinder according to Table 1.3 in Todreas & Kazimi.
Footnote (c): "Blends of fuel microspheres are molded to form fuel cylinders
each having a diameter of 15.7mm and a length of 5-6 cm."
For purposes of this homework, these cylinders will approximate an LWR
fuel pin with the fuel and clad radius both equal to 15.7 mm."""
htgr = ReactorType("HTGR", 2.3, None, df = 1.57, dco = 1.57,
                   shape = "hexagon", height = 630)
htgr.Q3 = 8.4
htgr.q1 = .0787
htgr.q1max = .298

# Advanced Gas Reactor
agr = ReactorType("AGR", ppitch = 2.57, npins = 37, df = 1.451, dco = 1.53,
                  shape = "circle", height = 829.6)
agr.Q3 = 2.66
agr.q1 = .17
agr.q1max = .298

# LMFBR
lmfbr = ReactorType("SFBR", ppitch = None, npins = None, df = 0.714, dco = .85,
                    shape = "hexagon", height = 270)
lmfbr.Q3 = 280
lmfbr.q1 = 0.29
lmfbr.q1max = 0.45
