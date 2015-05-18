from math import pi, cos

class elastanceController:

	def __init__(self):

		self.m_periodicTime = 0.0; #\todo think about this for restarts!
		self.m_timeToMaximumElastance = 0.2782;
		self.m_timeToRelax = 0.1391;
		self.m_minimumElastance = 4.10246e-3;
		self.m_maximumElastance = 3.0827e-1;
		self.m_heartPeriod = 0.86;


	def updateControl(self, currentParameterValue, delt):

		self.updatePeriodicTime(delt)
		elastance = self.getElastance(currentParameterValue)

		return elastance


	def updatePeriodicTime(self, delt):
	
		self.m_periodicTime = self.m_periodicTime + delt
		# Keep m_periodicTime in the range [0,m_heartPeriod):
		if self.m_periodicTime >= self.m_heartPeriod:
			self.m_periodicTime = self.m_periodicTime - self.m_heartPeriod

	def getElastance(self, currentParameterValue):
		# *** analytical elastance function from:
		#     pope, s. r.; ellwein, l. m.; zapata, c. l.; novak, v.; kelley, c. t. & olufsen, m. s.  
		#     estimation and identification of parameters in a lumped cerebrovascular model.
		#     math biosci eng, 2009, 6, 93-115

		# This is the elastance function. It's defined piecewise:
		if self.m_periodicTime <= self.m_timeToMaximumElastance:
			elastance = self.m_minimumElastance \
	        + 0.5*(self.m_maximumElastance - self.m_minimumElastance) \
	        * (1.0 - cos((self.m_periodicTime*pi)/self.m_timeToMaximumElastance))

		elif self.m_periodicTime <= (self.m_timeToMaximumElastance + self.m_timeToRelax):
		 	elastance = self.m_minimumElastance \
		    + 0.5*(self.m_maximumElastance-self.m_minimumElastance) \
		    * (1.0 + cos((self.m_periodicTime-self.m_timeToMaximumElastance)*(pi/self.m_timeToRelax)))

		elif self.m_periodicTime > (self.m_timeToMaximumElastance + self.m_timeToRelax):
			elastance = self.m_minimumElastance

		return elastance;
