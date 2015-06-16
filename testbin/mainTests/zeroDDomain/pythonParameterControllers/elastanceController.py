# def resolve_virtual_environment(override=None):
#     """Fetch the virtual environment path in the
#        process' environment or use an override."""
#     path = os.getenv('VIRTUAL_ENV')
#     if override:
#         path = os.path.join(os.getcwd(), override)
#     return path

# def activate_virtual_environment(environment_root):
#     """Configures the virtual environment starting at ``environment_root``."""
#     activate_script = os.path.join(
#         environment_root, 'Scripts', 'activate_this.py')
#     execfile(activate_script, {'__file__': activate_script})

# environment_root = resolve_virtual_environment(override)

# import sys
# sys.stderr = open('/home/carthurs/pyerrors.txt', 'w')
# # sys.path.pop()
# sys.path.append("/usr/lib/python2.7/")
# sys.path.append("/usr/lib/python2.7/dist-packages/")
# print sys.version_info
# sys.path.append('../')
# print sys.builtin_module_names
# import imp
# print imp.find_module('io')
# print imp.find_module('os')
# #print imp.find_module('numpy')
# # import io
# # sys.stdout = io.StringIO()
# # sys.stderr.write(sys.stdout.getvalue())
# # from math import *

# try:
# 	import numpy
# except ImportError as e:
# 	print "terrible thing happened!"
# 	print e.message
# 	print "terrible thing happened 2!"
# except:
# 	print "Unexpected error:", sys.exc_info()[0]

# try:
# 	import scipy
# except ImportError as e:
# 	print "terrible thing happened!"
# 	print e.message
# 	print "terrible thing happened 2!"
# except:
# 	print "Unexpected error:", sys.exc_info()[0]
# # import io
# # io = __import__('io', globals(), locals(), [], -1)
# # import logging
# sys.stderr.flush()
# # imp.load_module('scipy',None,'/usr/lib/python2.7/dist-packages/scipy', ('', '', 5))
# print sys.maxunicode
# print sys.path
#import io
#from random import getstate
#from numpy import *
#from scipy.integrate import ode
#import sys
#sys.path.append('/usr/lib64/python2.6/lib-dynload/')
from math import pi, cos

class elastanceController:

	def __init__(self):
		# import io

		self.m_periodicTime = 0.0; #\todo think about this for restarts!
		self.m_timeToMaximumElastance = 0.2782;
		self.m_timeToRelax = 0.1391;
		self.m_minimumElastance = 4.10246e-3;
		self.m_maximumElastance = 3.0827e-1;
		self.m_heartPeriod = 0.86;

	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByComponentIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):

		self.updatePeriodicTime(delt)
		elastance = self.getElastance(currentParameterValue)

		# for key in dictionaryOfPressuresByComponentIndex:
		# 	print "Pressure ", key, " was ", dictionaryOfPressuresByComponentIndex[key]
		# for key in dictionaryOfFlowsByComponentIndex:
		# 	print "Flow ", key, " was ", dictionaryOfFlowsByComponentIndex[key]
		# for key in dictionaryOfVolumesByComponentIndex:
		# 	print "Volume", key, "was", dictionaryOfVolumesByComponentIndex[key]

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
