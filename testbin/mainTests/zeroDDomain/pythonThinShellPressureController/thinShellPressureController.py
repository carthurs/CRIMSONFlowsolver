# Role of Individual Ionic Current Systems in Ventricular Cells Hypothesized by a Model Study, Matsuoka S, Sarai N, Kuratomi S, Ono K, and Noma A, 2003, The Japanese Journal of Physiology, 53, 105-123. PubMed ID: 12877767
# https://models.cellml.org/exposure/398d5dc7db9f2b9809abc29f440bd456/matsuoka_sarai_kuratomi_ono_noma_2003.cellml/view
import sys
sys.path.append("/usr/lib/python2.7/")
sys.path.append("/usr/lib/python2.7/dist-packages/")
from math import *
from numpy import *
from scipy.integrate import ode
from math import pi, cos
# try:
# 	import numpy
# except ImportError as e:
# 	print "terrible thing happened!"
# 	print e.message
# 	print "terrible thing happened 2!"
# except:
# 	print "Unexpected error:", sys.exc_info()[0]

class thinShellPressureController:

	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByComponentIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):

		self.time = self.time + delt * 1000
		# print "rem:",(self.time - self.stimStartTime) % self.stimPeriod
		# if ((self.time - self.stimStartTime) % self.stimPeriod == 0):
		# 	newPeriodJustStarting = True
		# else:
		# 	newPeriodJustStarting = False

		if self.time > 4000:
			self.stimPeriod = 600

		self.numberOfMuscleUnitsAroundCircumference = 205172 # Computed from Negroni96, equation (1), default constants[85] ( = L) and the circumference of a 100 ml sphere (to get Lm)
		circumference = self.computeCircumference(dictionaryOfVolumesByComponentIndex[5])
		self.halfSarcomereLength =  circumference / self.numberOfMuscleUnitsAroundCircumference * 1000 # *1000 is to convert to um # originally this was: constants[85] = 0.9623799975411884
		# self.halfSarcomereLength = halfSarcomereLength = 0.9623799975411884
		print "halfSarcomereLength:", self.halfSarcomereLength
		

		(init_states, constants) = self.initConsts(self.halfSarcomereLength)
		self.states = self.step_model(constants, self.states, circumference)

		self.computeAlgebraic(self.time, constants, self.states)
	# 	# elastance = max(self.algebraic[135] * self.m_maximumElastance, self.m_minimumElastance)

		# self.stepIndex = self.stepIndex + 1
		minimumTension = 0.0001
		tension = max(self.algebraic[135], minimumTension)
		print "tension:", tension
		print "self.algebraic[135]", self.algebraic[135]

		# print "elastance",elastance
		# print "time", self.time

		self.tensionHistory[self.stepIndex] = tension
		if self.time > 201:
			# reportedTension = (7*self.tensionHistory[self.stepIndex] +\
			# 					 2*self.tensionHistory[self.stepIndex-25] +\
			# 					 1*self.tensionHistory[self.stepIndex-50] +\
			# 					 2*self.tensionHistory[self.stepIndex-75] +\
			# 					 7*self.tensionHistory[self.stepIndex-100]) / 14.0
			# reportedTension = (8*self.tensionHistory[self.stepIndex] +\
			# 					 8*self.tensionHistory[self.stepIndex-25] +\
			# 					 8*self.tensionHistory[self.stepIndex-50] +\
			# 					 8*self.tensionHistory[self.stepIndex-75] +\
			# 					 8*self.tensionHistory[self.stepIndex-100]) / 40.0
			reportedTension = (8*self.tensionHistory[self.stepIndex] +\
								 8*self.tensionHistory[self.stepIndex-50] +\
								 8*self.tensionHistory[self.stepIndex-100] +\
								 8*self.tensionHistory[self.stepIndex-150] +\
								 8*self.tensionHistory[self.stepIndex-200]) / 40.0
			# reportedTension = (self.tensionHistory[self.stepIndex] +\
			# 					 self.tensionHistory[self.stepIndex-12] +\
			# 					 self.tensionHistory[self.stepIndex-25] +\
			# 					 self.tensionHistory[self.stepIndex-37] +\
			# 					 self.tensionHistory[self.stepIndex-50] +\
			# 					 self.tensionHistory[self.stepIndex-62] +\
			# 					 self.tensionHistory[self.stepIndex-75] +\
			# 					 self.tensionHistory[self.stepIndex-87] +\
			# 					 self.tensionHistory[self.stepIndex-100]) / 9.0
			# reportedTension = mean(self.tensionHistory[self.stepIndex-100:self.stepIndex+1])
			# reportedTension = (1*self.tensionHistory[self.stepIndex] +\
			# 					 2*self.tensionHistory[self.stepIndex-25] +\
			# 					 4*self.tensionHistory[self.stepIndex-50] +\
			# 					 6*self.tensionHistory[self.stepIndex-75] +\
			# 					 8*self.tensionHistory[self.stepIndex-100] +\
			# 					 6*self.tensionHistory[self.stepIndex-125] +\
			# 					 4*self.tensionHistory[self.stepIndex-150] +\
			# 					 2*self.tensionHistory[self.stepIndex-175] +\
			# 					 1*self.tensionHistory[self.stepIndex-200]) / 34.0
		else:
			reportedTension = minimumTension
		self.stepIndex = self.stepIndex + 1

		# reportedTension = tension
		print "reportedTension",reportedTension
		
		radius = self.computeRadius(dictionaryOfVolumesByComponentIndex[5]) # 100000
		# pressure = max(1000.0* self.computePressureFromLaplace(radius,reportedTension), 490.0)
		pressure = 1000.0* self.computePressureFromLaplace(radius,reportedTension)
		if pressure < 533:
			pressure = self.m_minimumElastance * dictionaryOfVolumesByComponentIndex[5]
			
		print "pressure returing:", pressure

		self.imposedPressureHistory[self.stepIndex] = pressure
		self.sarcomereHalfLengthHistory[self.stepIndex] = self.halfSarcomereLength
		self.tensionHistory[self.stepIndex] = tension
		self.customOutputArray1[self.stepIndex] = constants[85]# wobbles:self.algebraic[128]
		self.customOutputArray2[self.stepIndex] = self.states[36]# wobbles:self.algebraic[128]
		if self.time % 1000 == 0.0:
			savetxt('halfSarcHistory',self.sarcomereHalfLengthHistory)
			savetxt('pyPressHist',self.imposedPressureHistory)
			savetxt('tensionHistory',self.tensionHistory)
			savetxt('customOutput1',self.customOutputArray1)
			savetxt('customOutput2',self.customOutputArray2)
			


		# pressure = 533.3 * 20

		return pressure

	def computeCircumference(self,volume):
		radius = self.computeRadius(volume)
		circumference = 2 * pi * radius
		return circumference


	def computeRadius(self,volume):
		radius = power(array(volume * 3.0/4.0 / pi), 1/3.0)
		return radius

	def computePressureFromLaplace(self,radius,tension):
		pressure = 2.0 * tension * 8.0 / radius # 8.0 mm wall thickness
		return pressure


	# def updatePeriodicTime(self, delt):
	
	# 	self.m_periodicTime = self.m_periodicTime + delt
	# 	# Keep m_periodicTime in the range [0,m_heartPeriod):
	# 	if self.m_periodicTime >= self.m_heartPeriod:
	# 		self.m_periodicTime = self.m_periodicTime - self.m_heartPeriod

	# def getElastance(self, currentParameterValue):
	# 	# *** analytical elastance function from:
	# 	#     pope, s. self.r.; ellwein, l. m.; zapata, c. l.; novak, v.; kelley, c. t. & olufsen, m. s.  
	# 	#     estimation and identification of parameters in a lumped cerebrovascular model.
	# 	#     math biosci eng, 2009, 6, 93-115

	# 	# This is the elastance function. It's defined piecewise:
	# 	if self.m_periodicTime <= self.m_timeToMaximumElastance:
	# 		elastance = self.m_minimumElastance \
	#         + 0.5*(self.m_maximumElastance - self.m_minimumElastance) \
	#         * (1.0 - cos((self.m_periodicTime*pi)/self.m_timeToMaximumElastance))

	# 	elif self.m_periodicTime <= (self.m_timeToMaximumElastance + self.m_timeToRelax):
	# 	 	elastance = self.m_minimumElastance \
	# 	    + 0.5*(self.m_maximumElastance-self.m_minimumElastance) \
	# 	    * (1.0 + cos((self.m_periodicTime-self.m_timeToMaximumElastance)*(pi/self.m_timeToRelax)))

	# 	elif self.m_periodicTime > (self.m_timeToMaximumElastance + self.m_timeToRelax):
	# 		elastance = self.m_minimumElastance

	# 	return elastance;

	

	def createLegends(self):
	    legend_self.states = [""] * self.sizeself.states
	    legend_rates = [""] * self.sizeself.states
	    legend_algebraic = [""] * self.sizeAlgebraic
	    legend_voi = ""
	    legend_constants = [""] * self.sizeConstants
	    legend_voi = "time in component environment (millisecond)"
	    legend_self.states[0] = "Vm in component membrane (millivolt)"
	    legend_constants[0] = "self.r in component membrane (coulomb_millivolt_per_kelvin_millimole)"
	    legend_constants[1] = "T in component membrane (kelvin)"
	    legend_constants[2] = "F in component membrane (coulomb_per_millimole)"
	    legend_constants[3] = "Cm in component membrane (picoF)"
	    legend_algebraic[0] = "i_ext in component membrane (picoA)"
	    legend_algebraic[105] = "i_tot in component membrane (picoA)"
	    legend_algebraic[89] = "i_I in component membrane (picoA)"
	    legend_algebraic[45] = "i_Na in component sodium_current (picoA)"
	    legend_algebraic[50] = "i_Ca_L in component L_type_Ca_channel (picoA)"
	    legend_algebraic[54] = "i_Ca_T in component T_type_Ca_channel (picoA)"
	    legend_algebraic[70] = "i_K1 in component time_independent_potassium_current (picoA)"
	    legend_algebraic[71] = "i_Kr in component rapid_time_dependent_potassium_current (picoA)"
	    legend_algebraic[74] = "i_Ks in component slow_time_dependent_potassium_current (picoA)"
	    legend_algebraic[77] = "i_to in component transient_outward_current (picoA)"
	    legend_algebraic[103] = "i_NaK in component sodium_potassium_pump (picoA)"
	    legend_algebraic[94] = "i_NaCa in component sodium_calcium_exchanger (picoA)"
	    legend_algebraic[80] = "i_bNSC in component background_NSC_current (picoA)"
	    legend_algebraic[88] = "i_Cab in component background_Cab_current (picoA)"
	    legend_algebraic[81] = "i_Kpl in component background_Kpl_current (picoA)"
	    legend_algebraic[85] = "i_lCa in component background_lCa_current (picoA)"
	    legend_algebraic[87] = "i_KATP in component background_KATP_current (picoA)"
	    legend_constants[4] = "stim_start in component membrane (millisecond)"
	    legend_constants[5] = "stim_end in component membrane (millisecond)"
	    legend_constants[6] = "stim_period in component membrane (millisecond)"
	    legend_constants[7] = "stim_duration in component membrane (millisecond)"
	    legend_constants[8] = "stim_amplitude in component membrane (picoA)"
	    legend_constants[9] = "Nao in component external_ion_concentrations (millimolar)"
	    legend_constants[10] = "Cao in component external_ion_concentrations (millimolar)"
	    legend_constants[11] = "Ko in component external_ion_concentrations (millimolar)"
	    legend_self.states[1] = "Nai in component internal_ion_concentrations (millimolar)"
	    legend_algebraic[29] = "Cai in component internal_ion_concentrations (millimolar)"
	    legend_self.states[2] = "Ki in component internal_ion_concentrations (millimolar)"
	    legend_constants[12] = "Vi in component internal_ion_concentrations (micrometre3)"
	    legend_algebraic[106] = "i_net_Na in component internal_ion_concentrations (picoA)"
	    legend_algebraic[107] = "i_net_K in component internal_ion_concentrations (picoA)"
	    legend_algebraic[96] = "i_net_Ca in component internal_ion_concentrations (picoA)"
	    legend_algebraic[42] = "i_Na_Na in component sodium_current (picoA)"
	    legend_algebraic[48] = "i_CaL_Na in component L_type_Ca_channel (picoA)"
	    legend_algebraic[79] = "i_bNSC_Na in component background_NSC_current (picoA)"
	    legend_algebraic[84] = "i_lCa_Na in component background_lCa_current (picoA)"
	    legend_algebraic[75] = "i_to_K in component transient_outward_current (picoA)"
	    legend_algebraic[76] = "i_to_Na in component transient_outward_current (picoA)"
	    legend_algebraic[72] = "i_Ks_K in component slow_time_dependent_potassium_current (picoA)"
	    legend_algebraic[73] = "i_Ks_Na in component slow_time_dependent_potassium_current (picoA)"
	    legend_algebraic[44] = "i_Na_K in component sodium_current (picoA)"
	    legend_algebraic[49] = "i_CaL_K in component L_type_Ca_channel (picoA)"
	    legend_algebraic[78] = "i_bNSC_K in component background_NSC_current (picoA)"
	    legend_algebraic[83] = "i_lCa_K in component background_lCa_current (picoA)"
	    legend_algebraic[47] = "i_CaL_Ca in component L_type_Ca_channel (picoA)"
	    legend_algebraic[126] = "i_RyR in component RyR_channel (picoA)"
	    legend_algebraic[115] = "i_SR_U in component SR_calcium_pump (picoA)"
	    legend_algebraic[120] = "i_SR_L in component SR_L_current (picoA)"
	    legend_algebraic[140] = "dCaidt in component NL_model (millimolar_per_millisecond)"
	    legend_constants[13] = "CMDN_max in component internal_ion_concentrations (millimolar)"
	    legend_constants[14] = "K_mCMDN in component internal_ion_concentrations (millimolar)"
	    legend_self.states[3] = "Ca_Total in component internal_ion_concentrations (millimolar)"
	    legend_algebraic[13] = "b1 in component internal_ion_concentrations (millimolar)"
	    legend_algebraic[26] = "c1 in component internal_ion_concentrations (millimolar2)"
	    legend_algebraic[32] = "CF_Na in component constant_field_equations (millimolar)"
	    legend_algebraic[36] = "CF_Ca in component constant_field_equations (millimolar)"
	    legend_algebraic[39] = "CF_K in component constant_field_equations (millimolar)"
	    legend_self.states[4] = "ATPi in component ATP_production (millimolar)"
	    legend_algebraic[117] = "dATPdt in component NL_model (millimolar_per_millisecond)"
	    legend_constants[15] = "ProducingRate_Max in component ATP_production (per_millisecond)"
	    legend_constants[16] = "Adenosine_Total in component ATP_production (millimolar)"
	    legend_constants[17] = "P_Na in component sodium_current (picoA_per_millimolar)"
	    legend_self.states[5] = "p_AP_Na in component sodium_current_voltage_dependent_gate (dimensionless)"
	    legend_self.states[6] = "y in component sodium_current_ultra_slow_gate (dimensionless)"
	    legend_algebraic[1] = "p_RI_Na in component sodium_current_voltage_dependent_gate (dimensionless)"
	    legend_self.states[7] = "p_RP_Na in component sodium_current_voltage_dependent_gate (dimensionless)"
	    legend_self.states[8] = "p_AI_Na in component sodium_current_voltage_dependent_gate (dimensionless)"
	    legend_algebraic[14] = "k_RP_AP in component sodium_current_voltage_dependent_gate (per_millisecond)"
	    legend_algebraic[27] = "k_AP_RP in component sodium_current_voltage_dependent_gate (per_millisecond)"
	    legend_algebraic[37] = "k_RI_AI in component sodium_current_voltage_dependent_gate (per_millisecond)"
	    legend_algebraic[33] = "k_AI_RI in component sodium_current_voltage_dependent_gate (per_millisecond)"
	    legend_algebraic[30] = "k_AP_AI in component sodium_current_voltage_dependent_gate (per_millisecond)"
	    legend_constants[18] = "k_AI_AP in component sodium_current_voltage_dependent_gate (per_millisecond)"
	    legend_algebraic[40] = "k_RP_RI in component sodium_current_voltage_dependent_gate (per_millisecond)"
	    legend_algebraic[43] = "k_RI_RP in component sodium_current_voltage_dependent_gate (per_millisecond)"
	    legend_algebraic[2] = "alpha_y in component sodium_current_ultra_slow_gate (per_millisecond)"
	    legend_algebraic[15] = "beta_y in component sodium_current_ultra_slow_gate (per_millisecond)"
	    legend_algebraic[46] = "p_open_CaL in component L_type_Ca_channel (dimensionless)"
	    legend_algebraic[52] = "CaDiadic in component L_type_Ca_channel_Ca_dependent_gate (picoA)"
	    legend_constants[19] = "P_CaL in component L_type_Ca_channel (picoA_per_millimolar)"
	    legend_self.states[9] = "p_AP_CaL in component L_type_Ca_channel_voltage_dependent_gate (dimensionless)"
	    legend_self.states[10] = "p_U in component L_type_Ca_channel_Ca_dependent_gate (dimensionless)"
	    legend_self.states[11] = "p_UCa in component L_type_Ca_channel_Ca_dependent_gate (dimensionless)"
	    legend_self.states[12] = "y in component L_type_Ca_channel_ultra_slow_gate (dimensionless)"
	    legend_algebraic[3] = "p_RI_CaL in component L_type_Ca_channel_voltage_dependent_gate (dimensionless)"
	    legend_self.states[13] = "p_RP_CaL in component L_type_Ca_channel_voltage_dependent_gate (dimensionless)"
	    legend_self.states[14] = "p_AI_CaL in component L_type_Ca_channel_voltage_dependent_gate (dimensionless)"
	    legend_algebraic[16] = "k_RP_AP in component L_type_Ca_channel_voltage_dependent_gate (per_millisecond)"
	    legend_algebraic[28] = "k_AP_RP in component L_type_Ca_channel_voltage_dependent_gate (per_millisecond)"
	    legend_algebraic[31] = "k_RI_AI in component L_type_Ca_channel_voltage_dependent_gate (per_millisecond)"
	    legend_algebraic[34] = "k_AI_RI in component L_type_Ca_channel_voltage_dependent_gate (per_millisecond)"
	    legend_constants[20] = "k_AP_AI in component L_type_Ca_channel_voltage_dependent_gate (per_millisecond)"
	    legend_constants[21] = "k_AI_AP in component L_type_Ca_channel_voltage_dependent_gate (per_millisecond)"
	    legend_algebraic[38] = "k_RP_RI in component L_type_Ca_channel_voltage_dependent_gate (per_millisecond)"
	    legend_algebraic[41] = "k_RI_RP in component L_type_Ca_channel_voltage_dependent_gate (per_millisecond)"
	    legend_algebraic[51] = "iCaL in component L_type_Ca_channel_Ca_dependent_gate (picoA)"
	    legend_algebraic[53] = "Cacm in component L_type_Ca_channel_Ca_dependent_gate (millimolar)"
	    legend_algebraic[62] = "p_CCa in component L_type_Ca_channel_Ca_dependent_gate (dimensionless)"
	    legend_self.states[15] = "p_C in component L_type_Ca_channel_Ca_dependent_gate (dimensionless)"
	    legend_constants[22] = "k_CCa_UCa in component L_type_Ca_channel_Ca_dependent_gate (per_millisecond)"
	    legend_constants[23] = "k_UCa_CCa in component L_type_Ca_channel_Ca_dependent_gate (per_millisecond)"
	    legend_constants[24] = "k_C_U in component L_type_Ca_channel_Ca_dependent_gate (per_millisecond)"
	    legend_constants[25] = "k_U_C in component L_type_Ca_channel_Ca_dependent_gate (per_millisecond)"
	    legend_constants[92] = "k_UCa_U in component L_type_Ca_channel_Ca_dependent_gate (per_millisecond)"
	    legend_constants[26] = "k_U_UCa in component L_type_Ca_channel_Ca_dependent_gate (per_millimolar_millisecond)"
	    legend_constants[27] = "k_CCa_C in component L_type_Ca_channel_Ca_dependent_gate (per_millisecond)"
	    legend_constants[28] = "k_C_CCa in component L_type_Ca_channel_Ca_dependent_gate (per_millimolar_millisecond)"
	    legend_algebraic[55] = "CaEffC in component L_type_Ca_channel_Ca_dependent_gate (millimolar)"
	    legend_algebraic[57] = "CaEffU in component L_type_Ca_channel_Ca_dependent_gate (millimolar)"
	    legend_algebraic[60] = "k_UUCa_Ca in component L_type_Ca_channel_Ca_dependent_gate (per_millisecond)"
	    legend_algebraic[58] = "k_CCCa_Ca in component L_type_Ca_channel_Ca_dependent_gate (per_millisecond)"
	    legend_algebraic[4] = "alpha_y in component L_type_Ca_channel_ultra_slow_gate (per_millisecond)"
	    legend_algebraic[17] = "beta_y in component L_type_Ca_channel_ultra_slow_gate (per_millisecond)"
	    legend_constants[29] = "P_CaT in component T_type_Ca_channel (picoA_per_millimolar)"
	    legend_self.states[16] = "y1 in component T_type_Ca_channel_y1_gate (dimensionless)"
	    legend_self.states[17] = "y2 in component T_type_Ca_channel_y2_gate (dimensionless)"
	    legend_algebraic[5] = "alpha_y1 in component T_type_Ca_channel_y1_gate (per_millisecond)"
	    legend_algebraic[18] = "beta_y1 in component T_type_Ca_channel_y1_gate (per_millisecond)"
	    legend_algebraic[6] = "alpha_y2 in component T_type_Ca_channel_y2_gate (per_millisecond)"
	    legend_algebraic[19] = "beta_y2 in component T_type_Ca_channel_y2_gate (per_millisecond)"
	    legend_algebraic[56] = "E_K in component time_independent_potassium_current (millivolt)"
	    legend_constants[93] = "g_K1 in component time_independent_potassium_current (nanoS)"
	    legend_constants[30] = "P_K1_0 in component time_independent_potassium_current (nanoS_per_picoF)"
	    legend_algebraic[64] = "fO in component time_independent_potassium_current (dimensionless)"
	    legend_algebraic[65] = "fO2 in component time_independent_potassium_current (dimensionless)"
	    legend_algebraic[67] = "fO3 in component time_independent_potassium_current (dimensionless)"
	    legend_algebraic[69] = "fO4 in component time_independent_potassium_current (dimensionless)"
	    legend_algebraic[63] = "fB in component time_independent_potassium_current (dimensionless)"
	    legend_algebraic[59] = "mu in component time_independent_potassium_current (per_millisecond)"
	    legend_algebraic[61] = "lambda in component time_independent_potassium_current (per_millisecond)"
	    legend_self.states[18] = "y in component time_independent_potassium_current_y_gate (dimensionless)"
	    legend_algebraic[66] = "alpha_y in component time_independent_potassium_current_y_gate (per_millisecond)"
	    legend_algebraic[68] = "beta_y in component time_independent_potassium_current_y_gate (per_millisecond)"
	    legend_constants[94] = "g_Kr in component rapid_time_dependent_potassium_current (nanoS)"
	    legend_constants[31] = "P_Kr in component rapid_time_dependent_potassium_current (nanoS_per_picoF)"
	    legend_self.states[19] = "y1 in component rapid_time_dependent_potassium_current_y1_gate (dimensionless)"
	    legend_self.states[20] = "y2 in component rapid_time_dependent_potassium_current_y2_gate (dimensionless)"
	    legend_self.states[21] = "y3 in component rapid_time_dependent_potassium_current_y3_gate (dimensionless)"
	    legend_algebraic[7] = "alpha_y1 in component rapid_time_dependent_potassium_current_y1_gate (per_millisecond)"
	    legend_algebraic[20] = "beta_y1 in component rapid_time_dependent_potassium_current_y1_gate (per_millisecond)"
	    legend_algebraic[8] = "alpha_y2 in component rapid_time_dependent_potassium_current_y2_gate (per_millisecond)"
	    legend_algebraic[21] = "beta_y2 in component rapid_time_dependent_potassium_current_y2_gate (per_millisecond)"
	    legend_algebraic[9] = "alpha_y3 in component rapid_time_dependent_potassium_current_y3_gate (per_millisecond)"
	    legend_algebraic[22] = "beta_y3 in component rapid_time_dependent_potassium_current_y3_gate (per_millisecond)"
	    legend_self.states[22] = "y1 in component slow_time_dependent_potassium_current_y1_gate (dimensionless)"
	    legend_self.states[23] = "y2 in component slow_time_dependent_potassium_current_y2_gate (dimensionless)"
	    legend_constants[32] = "P_Ks_K in component slow_time_dependent_potassium_current (picoA_per_millimolar)"
	    legend_constants[33] = "P_Ks_Na in component slow_time_dependent_potassium_current (picoA_per_millimolar)"
	    legend_algebraic[10] = "alpha_y1 in component slow_time_dependent_potassium_current_y1_gate (per_millisecond)"
	    legend_algebraic[23] = "beta_y1 in component slow_time_dependent_potassium_current_y1_gate (per_millisecond)"
	    legend_algebraic[35] = "alpha_y2 in component slow_time_dependent_potassium_current_y2_gate (per_millisecond)"
	    legend_constants[34] = "beta_y2 in component slow_time_dependent_potassium_current_y2_gate (per_millisecond)"
	    legend_self.states[24] = "y1 in component transient_outward_current_y1_gate (dimensionless)"
	    legend_self.states[25] = "y2 in component transient_outward_current_y2_gate (dimensionless)"
	    legend_constants[35] = "P_to_K in component transient_outward_current (picoA_per_millimolar)"
	    legend_constants[36] = "P_to_Na in component transient_outward_current (picoA_per_millimolar)"
	    legend_algebraic[11] = "alpha_y1 in component transient_outward_current_y1_gate (per_millisecond)"
	    legend_algebraic[24] = "beta_y1 in component transient_outward_current_y1_gate (per_millisecond)"
	    legend_algebraic[12] = "alpha_y2 in component transient_outward_current_y2_gate (per_millisecond)"
	    legend_algebraic[25] = "beta_y2 in component transient_outward_current_y2_gate (per_millisecond)"
	    legend_constants[37] = "P_bNSC in component background_NSC_current (picoA_per_millimolar)"
	    legend_constants[95] = "P_Kpl in component background_Kpl_current (nanoS_per_millimolar)"
	    legend_constants[38] = "P_lCa in component background_lCa_current (picoA_per_millimolar)"
	    legend_algebraic[82] = "p_open in component background_lCa_current (dimensionless)"
	    legend_algebraic[86] = "p_open in component background_KATP_current (dimensionless)"
	    legend_constants[96] = "gamma in component background_KATP_current (nanoS)"
	    legend_constants[39] = "P_KATP in component background_KATP_current (nanoS_per_picoF)"
	    legend_constants[40] = "N in component background_KATP_current (picoF)"
	    legend_constants[41] = "P_Cab in component background_Cab_current (picoA_per_millimolar)"
	    legend_constants[97] = "p_E2Na in component sodium_calcium_exchanger (dimensionless)"
	    legend_algebraic[90] = "p_E1Na in component sodium_calcium_exchanger (dimensionless)"
	    legend_algebraic[91] = "p_E1Ca in component sodium_calcium_exchanger (dimensionless)"
	    legend_constants[100] = "p_E2Ca in component sodium_calcium_exchanger (dimensionless)"
	    legend_algebraic[92] = "k1 in component sodium_calcium_exchanger (per_millisecond)"
	    legend_algebraic[93] = "k2 in component sodium_calcium_exchanger (per_millisecond)"
	    legend_constants[42] = "k3 in component sodium_calcium_exchanger (per_millisecond)"
	    legend_constants[43] = "k4 in component sodium_calcium_exchanger (per_millisecond)"
	    legend_constants[44] = "Km_Nai in component sodium_calcium_exchanger (millimolar)"
	    legend_constants[45] = "Km_Nao in component sodium_calcium_exchanger (millimolar)"
	    legend_constants[46] = "Km_Cai in component sodium_calcium_exchanger (millimolar)"
	    legend_constants[47] = "Km_Cao in component sodium_calcium_exchanger (millimolar)"
	    legend_self.states[26] = "y in component sodium_calcium_exchanger_y_gate (dimensionless)"
	    legend_constants[48] = "P_NaCa in component sodium_calcium_exchanger (picoA_per_picoF)"
	    legend_constants[49] = "Partition in component sodium_calcium_exchanger (dimensionless)"
	    legend_algebraic[95] = "alpha_y in component sodium_calcium_exchanger_y_gate (per_millisecond)"
	    legend_algebraic[97] = "beta_y in component sodium_calcium_exchanger_y_gate (per_millisecond)"
	    legend_algebraic[102] = "p_E2Na in component sodium_potassium_pump (dimensionless)"
	    legend_algebraic[98] = "p_E1Na in component sodium_potassium_pump (dimensionless)"
	    legend_algebraic[99] = "p_E1K in component sodium_potassium_pump (dimensionless)"
	    legend_algebraic[104] = "p_E2K in component sodium_potassium_pump (dimensionless)"
	    legend_algebraic[100] = "k1 in component sodium_potassium_pump (per_millisecond)"
	    legend_constants[50] = "k2 in component sodium_potassium_pump (per_millisecond)"
	    legend_constants[51] = "k3 in component sodium_potassium_pump (per_millisecond)"
	    legend_constants[52] = "k4 in component sodium_potassium_pump (per_millisecond)"
	    legend_constants[53] = "Km_Nai in component sodium_potassium_pump (millimolar)"
	    legend_constants[54] = "Km_Nao in component sodium_potassium_pump (millimolar)"
	    legend_constants[55] = "Km_Ki in component sodium_potassium_pump (millimolar)"
	    legend_constants[56] = "Km_Ko in component sodium_potassium_pump (millimolar)"
	    legend_constants[57] = "Km_ATP in component sodium_potassium_pump (millimolar)"
	    legend_algebraic[101] = "Nao_Eff in component sodium_potassium_pump (millimolar)"
	    legend_self.states[27] = "y in component sodium_potassium_pump_y_gate (dimensionless)"
	    legend_constants[58] = "P_NaK in component sodium_potassium_pump (picoA_per_picoF)"
	    legend_algebraic[109] = "alpha_y in component sodium_potassium_pump_y_gate (per_millisecond)"
	    legend_algebraic[111] = "beta_y in component sodium_potassium_pump_y_gate (per_millisecond)"
	    legend_algebraic[110] = "p_E2Ca in component SR_calcium_pump (dimensionless)"
	    legend_algebraic[108] = "p_E1Ca in component SR_calcium_pump (dimensionless)"
	    legend_algebraic[112] = "p_E1 in component SR_calcium_pump (dimensionless)"
	    legend_algebraic[113] = "p_E2 in component SR_calcium_pump (dimensionless)"
	    legend_constants[59] = "k1 in component SR_calcium_pump (per_millisecond)"
	    legend_algebraic[114] = "k2 in component SR_calcium_pump (per_millisecond)"
	    legend_constants[60] = "k3 in component SR_calcium_pump (per_millisecond)"
	    legend_constants[61] = "k4 in component SR_calcium_pump (per_millisecond)"
	    legend_constants[62] = "Km_CaSR in component SR_calcium_pump (millimolar)"
	    legend_constants[63] = "Km_CaCyto in component SR_calcium_pump (millimolar)"
	    legend_constants[64] = "Km_ATP in component SR_calcium_pump (millimolar)"
	    legend_constants[65] = "i_max in component SR_calcium_pump (picoA)"
	    legend_self.states[28] = "Caup in component Ca_concentrations_in_SR (millimolar)"
	    legend_self.states[29] = "y in component SR_calcium_pump_y_gate (dimensionless)"
	    legend_algebraic[116] = "alpha_y in component SR_calcium_pump_y_gate (per_millisecond)"
	    legend_algebraic[118] = "beta_y in component SR_calcium_pump_y_gate (per_millisecond)"
	    legend_constants[66] = "P_RyR in component RyR_channel (picoA_per_millimolar)"
	    legend_algebraic[119] = "k1 in component RyR_channel (per_millisecond)"
	    legend_algebraic[124] = "k2 in component RyR_channel (per_millisecond)"
	    legend_algebraic[125] = "k3 in component RyR_channel (per_millisecond)"
	    legend_constants[67] = "k4 in component RyR_channel (per_millisecond)"
	    legend_self.states[30] = "p_open_RyR in component RyR_channel (dimensionless)"
	    legend_self.states[31] = "p_close_RyR in component RyR_channel (dimensionless)"
	    legend_algebraic[123] = "Carel in component Ca_concentrations_in_SR (millimolar)"
	    legend_constants[68] = "Diadid_Factor in component RyR_channel (per_picoA_millisecond)"
	    legend_algebraic[127] = "i_SR_T in component SR_T_current (picoA)"
	    legend_constants[69] = "P_SR_T in component SR_T_current (picoA_per_millimolar)"
	    legend_constants[70] = "P_SR_L in component SR_L_current (picoA_per_millimolar)"
	    legend_self.states[32] = "Ca_Total in component Ca_concentrations_in_SR (millimolar)"
	    legend_constants[71] = "V_rel in component Ca_concentrations_in_SR (micrometre3)"
	    legend_constants[72] = "V_up in component Ca_concentrations_in_SR (micrometre3)"
	    legend_constants[73] = "CSQN_max in component Ca_concentrations_in_SR (millimolar)"
	    legend_constants[74] = "K_mCSQN in component Ca_concentrations_in_SR (millimolar)"
	    legend_algebraic[121] = "b1 in component Ca_concentrations_in_SR (millimolar)"
	    legend_algebraic[122] = "c1 in component Ca_concentrations_in_SR (millimolar2)"
	    legend_constants[98] = "EffFraction in component NL_model (dimensionless)"
	    legend_self.states[33] = "pCa in component NL_model (dimensionless)"
	    legend_self.states[34] = "pCaCB in component NL_model (dimensionless)"
	    legend_self.states[35] = "pCB in component NL_model (dimensionless)"
	    legend_algebraic[130] = "p in component NL_model (dimensionless)"
	    legend_constants[75] = "T_t in component NL_model (millimolar)"
	    legend_algebraic[134] = "Q_a in component NL_model (per_millisecond)"
	    legend_algebraic[132] = "Q_b in component NL_model (per_millisecond)"
	    legend_algebraic[136] = "Q_r in component NL_model (per_millisecond)"
	    legend_algebraic[137] = "Q_d in component NL_model (per_millisecond)"
	    legend_algebraic[138] = "Q_d1 in component NL_model (per_millisecond)"
	    legend_algebraic[139] = "Q_d2 in component NL_model (per_millisecond)"
	    legend_constants[76] = "Y_1 in component NL_model (per_millimolar_millisecond)"
	    legend_constants[77] = "Y_2 in component NL_model (per_millisecond)"
	    legend_constants[78] = "Y_3 in component NL_model (per_millisecond)"
	    legend_constants[79] = "Y_4 in component NL_model (per_millisecond)"
	    legend_constants[80] = "Y_d in component NL_model (millisecond_per_micrometre2)"
	    legend_constants[81] = "Z_1 in component NL_model (per_millisecond)"
	    legend_constants[82] = "Z_2 in component NL_model (per_millisecond)"
	    legend_constants[83] = "Z_3 in component NL_model (per_millimolar_millisecond)"
	    legend_algebraic[128] = "h in component NL_model (micrometre)"
	    legend_constants[84] = "L_a in component NL_model (micrometre)"
	    legend_constants[85] = "L in component NL_model (micrometre)"
	    legend_algebraic[133] = "ForceCB in component NL_model (mN_per_mm2)"
	    legend_self.states[36] = "X in component NL_model (micrometre)"
	    legend_algebraic[131] = "NewCBF in component NL_model (mN_per_mm2_micrometre)"
	    legend_algebraic[129] = "CBBound in component NL_model (millimolar)"
	    legend_constants[86] = "KForceEC in component NL_model (mN_per_mm2_micrometre5)"
	    legend_constants[87] = "ZeroForceEL in component NL_model (micrometre)"
	    legend_constants[88] = "KForceLinearEc in component NL_model (mN_per_mm2_micrometre)"
	    legend_constants[89] = "ForceFactor in component NL_model (mN_per_mm2_micrometre_millimolar)"
	    legend_constants[99] = "ForceEcomp in component NL_model (mN_per_mm2)"
	    legend_constants[90] = "B in component NL_model (per_millisecond)"
	    legend_constants[91] = "h_c in component NL_model (micrometre)"
	    legend_algebraic[135] = "ForceExt in component NL_model (mN_per_mm2)"
	    legend_rates[0] = "d/dt Vm in component membrane (millivolt)"
	    legend_rates[1] = "d/dt Nai in component internal_ion_concentrations (millimolar)"
	    legend_rates[2] = "d/dt Ki in component internal_ion_concentrations (millimolar)"
	    legend_rates[3] = "d/dt Ca_Total in component internal_ion_concentrations (millimolar)"
	    legend_rates[4] = "d/dt ATPi in component ATP_production (millimolar)"
	    legend_rates[7] = "d/dt p_RP_Na in component sodium_current_voltage_dependent_gate (dimensionless)"
	    legend_rates[5] = "d/dt p_AP_Na in component sodium_current_voltage_dependent_gate (dimensionless)"
	    legend_rates[8] = "d/dt p_AI_Na in component sodium_current_voltage_dependent_gate (dimensionless)"
	    legend_rates[6] = "d/dt y in component sodium_current_ultra_slow_gate (dimensionless)"
	    legend_rates[13] = "d/dt p_RP_CaL in component L_type_Ca_channel_voltage_dependent_gate (dimensionless)"
	    legend_rates[9] = "d/dt p_AP_CaL in component L_type_Ca_channel_voltage_dependent_gate (dimensionless)"
	    legend_rates[14] = "d/dt p_AI_CaL in component L_type_Ca_channel_voltage_dependent_gate (dimensionless)"
	    legend_rates[10] = "d/dt p_U in component L_type_Ca_channel_Ca_dependent_gate (dimensionless)"
	    legend_rates[11] = "d/dt p_UCa in component L_type_Ca_channel_Ca_dependent_gate (dimensionless)"
	    legend_rates[15] = "d/dt p_C in component L_type_Ca_channel_Ca_dependent_gate (dimensionless)"
	    legend_rates[12] = "d/dt y in component L_type_Ca_channel_ultra_slow_gate (dimensionless)"
	    legend_rates[16] = "d/dt y1 in component T_type_Ca_channel_y1_gate (dimensionless)"
	    legend_rates[17] = "d/dt y2 in component T_type_Ca_channel_y2_gate (dimensionless)"
	    legend_rates[18] = "d/dt y in component time_independent_potassium_current_y_gate (dimensionless)"
	    legend_rates[19] = "d/dt y1 in component rapid_time_dependent_potassium_current_y1_gate (dimensionless)"
	    legend_rates[20] = "d/dt y2 in component rapid_time_dependent_potassium_current_y2_gate (dimensionless)"
	    legend_rates[21] = "d/dt y3 in component rapid_time_dependent_potassium_current_y3_gate (dimensionless)"
	    legend_rates[22] = "d/dt y1 in component slow_time_dependent_potassium_current_y1_gate (dimensionless)"
	    legend_rates[23] = "d/dt y2 in component slow_time_dependent_potassium_current_y2_gate (dimensionless)"
	    legend_rates[24] = "d/dt y1 in component transient_outward_current_y1_gate (dimensionless)"
	    legend_rates[25] = "d/dt y2 in component transient_outward_current_y2_gate (dimensionless)"
	    legend_rates[26] = "d/dt y in component sodium_calcium_exchanger_y_gate (dimensionless)"
	    legend_rates[27] = "d/dt y in component sodium_potassium_pump_y_gate (dimensionless)"
	    legend_rates[29] = "d/dt y in component SR_calcium_pump_y_gate (dimensionless)"
	    legend_rates[30] = "d/dt p_open_RyR in component RyR_channel (dimensionless)"
	    legend_rates[31] = "d/dt p_close_RyR in component RyR_channel (dimensionless)"
	    legend_rates[32] = "d/dt Ca_Total in component Ca_concentrations_in_SR (millimolar)"
	    legend_rates[28] = "d/dt Caup in component Ca_concentrations_in_SR (millimolar)"
	    legend_rates[36] = "d/dt X in component NL_model (micrometre)"
	    legend_rates[33] = "d/dt pCa in component NL_model (dimensionless)"
	    legend_rates[34] = "d/dt pCaCB in component NL_model (dimensionless)"
	    legend_rates[35] = "d/dt pCB in component NL_model (dimensionless)"
	    return (legend_self.states, legend_algebraic, legend_voi, legend_constants)

	def initConsts(self,halfSarcomereLength):
		constants = array([0.0] * self.sizeConstants); states = [0.0] * self.sizeStates;
		# print "self.sizeStates:", self.sizeStates
		states[0] = -85.95752434460744
		constants[0] = 8.3143
		constants[1] = 310
		constants[2] = 96.4867
		constants[3] = 132
		constants[4] = self.stimStartTime # 50
		constants[5] = 1000000
		constants[6] = self.stimPeriod # 400
		constants[7] = self.stimDuration #2
		constants[8] = -4000
		constants[9] = 140
		constants[10] = 1.8
		constants[11] = 5.4
		states[1] = 4.925761439682025
		states[2] = 143.1837333000449
		constants[12] = 8000
		constants[13] = 0.05
		constants[14] = 0.00238
		states[3] = 4.0180173572968586e-4
		states[4] = 4.657102729020499
		constants[15] = 0.003
		constants[16] = 5
		constants[17] = 2860
		states[5] = 1.779648367445368e-5
		states[6] = 0.5861887862983165
		states[7] = 0.3556412697995689
		states[8] = 0.40285968661346977
		constants[18] = 0.0000875
		constants[19] = 8712
		states[9] = 1.5445004166497696e-6
		states[10] = 0.17246483915629204
		states[11] = 6.098246017787626e-5
		states[12] = 0.9985266538252986
		states[13] = 0.9968480629364956
		states[14] = 8.77325391245903e-4
		constants[20] = 0.004
		constants[21] = 0.001
		states[15] = 0.4250747299372254
		constants[22] = 0.0003
		constants[23] = 0.35
		constants[24] = 0.143
		constants[25] = 0.35
		constants[26] = 6.954
		constants[27] = 0.0042
		constants[28] = 6.954
		constants[29] = 612
		states[16] = 1.6882718240109127e-5
		states[17] = 0.8585352091865849
		constants[30] = 1.146
		states[18] = 0.6080573900752752
		constants[31] = 0.00864
		states[19] = 0.0018339931180983765
		states[20] = 0.20443083454225305
		states[21] = 0.967887666264921
		states[22] = 0.09738789658609195
		states[23] = 0.09745345578743213
		constants[32] = 5.04
		constants[33] = 0.2016
		constants[34] = 0.004444
		states[24] = 7.956883250874798e-4
		states[25] = 0.9999125083105881
		constants[35] = 0.033
		constants[36] = 0.00297
		constants[37] = 0.385
		constants[38] = 0.11
		constants[39] = 0.0236
		constants[40] = 2333
		constants[41] = 0.04
		constants[42] = 1
		constants[43] = 1
		constants[44] = 8.75
		constants[45] = 87.5
		constants[46] = 0.00138
		constants[47] = 1.38
		states[26] = 0.9891789193465331
		constants[48] = 6.81
		constants[49] = 0.32
		constants[50] = 0.04
		constants[51] = 0.01
		constants[52] = 0.165
		constants[53] = 4.05
		constants[54] = 69.8
		constants[55] = 32.88
		constants[56] = 0.258
		constants[57] = 0.094
		states[27] = 0.5910747147428818
		constants[58] = 21
		constants[59] = 0.01
		constants[60] = 1
		constants[61] = 0.01
		constants[62] = 0.08
		constants[63] = 0.0008
		constants[64] = 0.1
		constants[65] = 162500
		states[28] = 2.611712901567567
		states[29] = 0.46108441538480216
		constants[66] = 62000
		constants[67] = 0.000849
		states[30] = 3.4314360001543243e-4
		states[31] = 0.19135178123107768
		constants[68] = -150
		constants[69] = 386
		constants[70] = 459
		states[32] = 9.455741736977666
		constants[71] = 160
		constants[72] = 400
		constants[73] = 10
		constants[74] = 0.8
		states[33] = 0.02490898775497523
		states[34] = 0.001990153835322864
		states[35] = 4.2941813853474524e-4
		constants[75] = 0.07
		constants[76] = 39
		constants[77] = 0.0039
		constants[78] = 0.03
		constants[79] = 0.12
		constants[80] = 0.027
		constants[81] = 0.03
		constants[82] = 0.0039
		constants[83] = 1560
		constants[84] = 1.17
		constants[85] = halfSarcomereLength #0.9623799975411884
		states[36] = 0.9573749975411884
		constants[86] = 140000
		constants[87] = 0.97
		constants[88] = 200
		constants[89] = 1800000
		constants[90] = 1.2
		constants[91] = 0.005
		constants[92] = (constants[27]*constants[	24]*constants[26]*constants[23])/(constants[25]*constants[28]*constants[22])
		constants[93] = constants[30]*constants[3]*(power(constants[11]/5.40000, 0.400000))
		constants[94] = constants[31]*constants[3]*(power(constants[11]/5.40000, 0.200000))
		constants[95] = 0.000110000*(power(constants[11]/5.40000, 0.160000))
		constants[96] = constants[39]*constants[40]*(power(constants[11]/1.00000, 0.240000))
		constants[97] = 1.00000/(1.00000+(power(constants[45]/constants[9], 3.00000))*(1.00000+constants[10]/constants[47]))
		constants[98] = exp(-20.0000*(power(constants[85]-constants[84], 2.00000)))
		# constants[99] = constants[86]*(power(constants[87]-constants[85], 5.00000))+constants[88]*(constants[87]-constants[85])
		constants[99] = select([constants[87]-constants[85] >= 0.0, constants[87]-constants[85] < 0.0],[30.0*(power(constants[87]-constants[85], 1.00000))+constants[88]*(constants[87]-constants[85]), constants[86]*(power(constants[87]-constants[85], 5.00000))+constants[88]*(constants[87]-constants[85])])
		constants[100] = 1.00000/(1.00000+(constants[47]/constants[10])*(1.00000+power(constants[9]/constants[45], 3.00000)))
		return (states, constants)

	def computeRates(self, voi, states, constants):
	    self.rates = [0.0] * self.sizeStates; self.algebraic = [0.0] * self.sizeAlgebraic
	    self.algebraic[2] = 1.00000/(9.00000e+09*exp(states[0]/5.00000)+8000.00*exp(states[0]/100.000))
	    self.algebraic[15] = 1.00000/(0.0140000*exp(-states[0]/5.00000)+4000.00*exp(-states[0]/100.000))
	    self.rates[6] = self.algebraic[2]*(1.00000-states[6])-self.algebraic[15]*states[6]
	    self.algebraic[4] = 1.00000/(250000.*exp(states[0]/9.00000)+58.0000*exp(states[0]/65.0000))
	    self.algebraic[17] = 1.00000/(1800.00*exp(-states[0]/14.0000)+66.0000*exp(-states[0]/65.0000))
	    self.rates[12] = self.algebraic[4]*(1.00000-states[12])-self.algebraic[17]*states[12]
	    self.algebraic[5] = 1.00000/(0.0190000*exp(-states[0]/5.60000)+0.820000*exp(-states[0]/250.000))
	    self.algebraic[18] = 1.00000/(40.0000*exp(states[0]/6.30000)+1.50000*exp(states[0]/10000.0))
	    self.rates[16] = self.algebraic[5]*(1.00000-states[16])-self.algebraic[18]*states[16]
	    self.algebraic[6] = 1.00000/(62000.0*exp(states[0]/10.1000)+30.0000*exp(states[0]/3000.00))
	    self.algebraic[19] = 1.00000/(0.000600000*exp(-states[0]/6.70000)+1.20000*exp(-states[0]/25.0000))
	    self.rates[17] = self.algebraic[6]*(1.00000-states[17])-self.algebraic[19]*states[17]
	    self.algebraic[7] = 1.00000/(20.0000*exp(-states[0]/11.5000)+5.00000*exp(-states[0]/300.000))
	    self.algebraic[20] = 1.00000/(160.000*exp(states[0]/28.0000)+200.000*exp(states[0]/1000.00))+1.00000/(2500.00*exp(states[0]/20.0000))
	    self.rates[19] = self.algebraic[7]*(1.00000-states[19])-self.algebraic[20]*states[19]
	    self.algebraic[8] = 1.00000/(200.000*exp(-states[0]/13.0000)+20.0000*exp(-states[0]/300.000))
	    self.algebraic[21] = 1.00000/(1600.00*exp(states[0]/28.0000)+2000.00*exp(states[0]/1000.00))+1.00000/(10000.0*exp(states[0]/20.0000))
	    self.rates[20] = self.algebraic[8]*(1.00000-states[20])-self.algebraic[21]*states[20]
	    self.algebraic[9] = 1.00000/(10.0000*exp(states[0]/17.0000)+2.50000*exp(states[0]/300.000))
	    self.algebraic[22] = 1.00000/(0.350000*exp(-states[0]/17.0000)+2.00000*exp(-states[0]/150.000))
	    self.rates[21] = self.algebraic[9]*(1.00000-states[21])-self.algebraic[22]*states[21]
	    self.algebraic[10] = 1.00000/(85.0000*exp(-states[0]/10.5000)+370.000*exp(-states[0]/62.0000))
	    self.algebraic[23] = 1.00000/(1450.00*exp(states[0]/20.0000)+260.000*exp(states[0]/100.000))
	    self.rates[22] = self.algebraic[10]*(1.00000-states[22])-self.algebraic[23]*states[22]
	    self.algebraic[11] = 1.00000/(11.0000*exp(-states[0]/28.0000)+0.200000*exp(-states[0]/400.000))
	    self.algebraic[24] = 1.00000/(4.40000*exp(states[0]/16.0000)+0.200000*exp(states[0]/500.000))
	    self.rates[24] = self.algebraic[11]*(1.00000-states[24])-self.algebraic[24]*states[24]
	    self.algebraic[12] = (0.00380000*exp(-(states[0]+13.5000)/11.3000))/(1.00000+0.0513350*exp(-(states[0]+13.5000)/11.3000))
	    self.algebraic[25] = (0.00380000*exp((states[0]+13.5000)/11.3000))/(1.00000+0.0670830*exp((states[0]+13.5000)/11.3000))
	    self.rates[25] = self.algebraic[12]*(1.00000-states[25])-self.algebraic[25]*states[25]
	    self.algebraic[16] = 1.00000/(0.270000*exp(-states[0]/5.90000)+1.50000*exp(-states[0]/65.0000))
	    self.algebraic[28] = 1.00000/(480.000*exp(states[0]/7.00000)+2.20000*exp(states[0]/65.0000))
	    self.rates[9] = (states[13]*self.algebraic[16]+states[14]*constants[21])-states[9]*(self.algebraic[28]+constants[20])
	    self.algebraic[14] = 1.00000/(0.102700*exp(-states[0]/8.00000)+0.250000*exp(-states[0]/50.0000))
	    self.algebraic[27] = 1.00000/(26.0000*exp(states[0]/17.0000)+0.0200000*exp(states[0]/800.000))
	    self.algebraic[30] = 1.00000/(0.800000*exp(-states[0]/400.000))
	    self.rates[5] = (states[7]*self.algebraic[14]+states[8]*constants[18])-states[5]*(self.algebraic[27]+self.algebraic[30])
	    self.algebraic[3] = ((1.00000-states[9])-states[13])-states[14]
	    self.algebraic[31] = 1.00000/(0.00180000*exp(-states[0]/7.40000)+2.00000*exp(-states[0]/100.000))
	    self.algebraic[34] = 1.00000/(2.20000e+06*exp(states[0]/7.40000)+11.0000*exp(states[0]/100.000))
	    self.rates[14] = (self.algebraic[3]*self.algebraic[31]+states[9]*constants[20])-states[14]*(self.algebraic[34]+constants[21])
	    self.algebraic[13] = (constants[13]-states[3])+constants[14]
	    self.algebraic[26] = constants[14]*states[3]
	    self.algebraic[29] = (power(power(self.algebraic[13], 2.00000)+4.00000*self.algebraic[26], 1.0/2)-self.algebraic[13])/2.00000
	    self.algebraic[35] = 3.70000*self.algebraic[29]
	    self.rates[23] = self.algebraic[35]*(1.00000-states[23])-constants[34]*states[23]
	    self.algebraic[1] = ((1.00000-states[7])-states[5])-states[8]
	    self.algebraic[37] = 1.00000/(0.000102700*exp(-states[0]/8.00000)+5.00000*exp(-states[0]/400.000))
	    self.algebraic[33] = 1.00000/(1300.00*exp(states[0]/20.0000)+0.0400000*exp(states[0]/800.000))
	    self.rates[8] = (self.algebraic[1]*self.algebraic[37]+states[5]*self.algebraic[30])-states[8]*(self.algebraic[33]+constants[18])
	    self.algebraic[38] = 0.0400000/(1.00000+(constants[21]*self.algebraic[28]*self.algebraic[31])/(constants[20]*self.algebraic[16]*self.algebraic[34]))
	    self.algebraic[41] = 0.0400000-self.algebraic[38]
	    self.rates[13] = (states[9]*self.algebraic[28]+self.algebraic[3]*self.algebraic[41])-states[13]*(self.algebraic[38]+self.algebraic[16])
	    self.algebraic[40] = 0.0100000/(1.00000+(constants[18]*self.algebraic[27]*self.algebraic[37])/(self.algebraic[30]*self.algebraic[14]*self.algebraic[33]))
	    self.algebraic[43] = 0.0100000-self.algebraic[40]
	    self.rates[7] = (states[5]*self.algebraic[27]+self.algebraic[1]*self.algebraic[43])-states[7]*(self.algebraic[40]+self.algebraic[14])
	    self.algebraic[36] = self.custom_piecewise([equal(states[0] , 0.00000), -constants[10] , True, (((2.00000*constants[2]*states[0])/(constants[0]*constants[1]))*(self.algebraic[29]-constants[10]*exp((-2.00000*constants[2]*states[0])/(constants[0]*constants[1]))))/(1.00000-exp((-2.00000*constants[2]*states[0])/(constants[0]*constants[1])))])
	    self.algebraic[51] = 0.0676000*self.algebraic[36]
	    self.algebraic[53] = self.algebraic[29]-0.300000*self.algebraic[51]
	    self.algebraic[55] = self.algebraic[53]*states[9]
	    self.algebraic[57] = self.algebraic[55]+self.algebraic[29]*(1.00000-states[9])
	    self.algebraic[60] = constants[26]*self.algebraic[57]
	    self.rates[10] = (states[15]*constants[24]+states[11]*constants[92])-states[10]*(self.algebraic[60]+constants[25])
	    self.algebraic[62] = ((1.00000-states[15])-states[10])-states[11]
	    self.rates[11] = (states[10]*self.algebraic[60]+self.algebraic[62]*constants[22])-states[11]*(constants[23]+constants[92])
	    self.rates[15] = (self.algebraic[62]*constants[27]+states[10]*constants[25])-states[15]*(constants[24]+constants[28]*self.algebraic[53]*states[9])
	    self.algebraic[56] = ((constants[0]*constants[1])/constants[2])*log(constants[11]/states[2])
	    self.algebraic[66] = 1.00000/(8000.00*exp(((states[0]-self.algebraic[56])-97.0000)/8.50000)+7.00000*exp(((states[0]-self.algebraic[56])-97.0000)/300.000))
	    self.algebraic[59] = (0.750000*exp(0.0350000*((states[0]-self.algebraic[56])-10.0000)))/(1.00000+exp(0.0150000*((states[0]-self.algebraic[56])-140.000)))
	    self.algebraic[61] = (3.00000*exp(-0.0480000*((states[0]-self.algebraic[56])-10.0000))*(1.00000+exp(0.0640000*((states[0]-self.algebraic[56])-38.0000))))/(1.00000+exp(0.0300000*((states[0]-self.algebraic[56])-70.0000)))
	    self.algebraic[64] = self.algebraic[61]/(self.algebraic[59]+self.algebraic[61])
	    self.algebraic[68] = ((power(self.algebraic[64], 4.00000))*1.00000)/(0.000140000*exp(-((states[0]-self.algebraic[56])-97.0000)/9.10000)+0.200000*exp(-((states[0]-self.algebraic[56])-97.0000)/500.000))
	    self.rates[18] = self.algebraic[66]*(1.00000-states[18])-self.algebraic[68]*states[18]
	    self.algebraic[93] = 1.00000*exp(((constants[49]-1.00000)*constants[2]*states[0])/(constants[0]*constants[1]))
	    self.algebraic[95] = self.algebraic[93]*constants[97]+constants[43]*constants[100]
	    self.algebraic[90] = 1.00000/(1.00000+(power(constants[44]/states[1], 3.00000))*(1.00000+self.algebraic[29]/constants[46]))
	    self.algebraic[91] = 1.00000/(1.00000+(constants[46]/self.algebraic[29])*(1.00000+power(states[1]/constants[44], 3.00000)))
	    self.algebraic[92] = 1.00000*exp((constants[49]*constants[2]*states[0])/(constants[0]*constants[1]))
	    self.algebraic[97] = self.algebraic[92]*self.algebraic[90]+constants[42]*self.algebraic[91]
	    self.rates[26] = self.algebraic[95]*(1.00000-states[26])-self.algebraic[97]*states[26]
	    self.algebraic[0] = self.custom_piecewise([greater_equal(voi , constants[4]) & less_equal(voi , constants[5]) & less_equal((voi-constants[4])-floor((voi-constants[4])/constants[6])*constants[6] , constants[7]), constants[8] , True, 0.00000])
	    self.algebraic[32] = self.custom_piecewise([equal(states[0] , 0.00000), -constants[9] , True, (((constants[2]*states[0])/(constants[0]*constants[1]))*(states[1]-constants[9]*exp((-constants[2]*states[0])/(constants[0]*constants[1]))))/(1.00000-exp((-constants[2]*states[0])/(constants[0]*constants[1])))])
	    self.algebraic[79] = constants[37]*self.algebraic[32]
	    self.algebraic[39] = self.custom_piecewise([equal(states[0] , 0.00000), states[2] , True, (((constants[2]*states[0])/(constants[0]*constants[1]))*(states[2]-constants[11]*exp((-constants[2]*states[0])/(constants[0]*constants[1]))))/(1.00000-exp((-constants[2]*states[0])/(constants[0]*constants[1])))])
	    self.algebraic[78] = 0.400000*constants[37]*self.algebraic[39]
	    self.algebraic[80] = self.algebraic[78]+self.algebraic[79]
	    self.algebraic[88] = constants[41]*self.algebraic[36]
	    self.algebraic[81] = self.custom_piecewise([equal(states[0] , -3.00000), constants[95]*self.algebraic[39]*13.0077 , True, (constants[95]*self.algebraic[39]*(states[0]+3.00000))/(1.00000-exp(-(states[0]+3.00000)/13.0000))])
	    self.algebraic[82] = 1.00000/(1.00000+power(0.00120000/self.algebraic[29], 3.00000))
	    self.algebraic[84] = constants[38]*self.algebraic[32]*self.algebraic[82]
	    self.algebraic[83] = constants[38]*self.algebraic[39]*self.algebraic[82]
	    self.algebraic[85] = self.algebraic[83]+self.algebraic[84]
	    self.algebraic[86] = 0.800000/(1.00000+power(states[4]/0.100000, 2.00000))
	    self.algebraic[87] = constants[96]*(states[0]-self.algebraic[56])*self.algebraic[86]
	    self.algebraic[89] = self.algebraic[80]+self.algebraic[88]+self.algebraic[81]+self.algebraic[85]+self.algebraic[87]
	    self.algebraic[42] = constants[17]*self.algebraic[32]*states[5]*states[6]
	    self.algebraic[44] = 0.100000*constants[17]*self.algebraic[39]*states[5]*states[6]
	    self.algebraic[45] = self.algebraic[42]+self.algebraic[44]
	    self.algebraic[46] = (states[9]*(states[10]+states[11])*states[12])/(1.00000+power(1.40000/states[4], 3.00000))
	    self.algebraic[48] = 1.85000e-05*constants[19]*self.algebraic[32]*self.algebraic[46]
	    self.algebraic[49] = 0.000365000*constants[19]*self.algebraic[39]*self.algebraic[46]
	    self.algebraic[47] = constants[19]*self.algebraic[36]*self.algebraic[46]
	    self.algebraic[50] = self.algebraic[48]+self.algebraic[47]+self.algebraic[49]
	    self.algebraic[54] = constants[29]*self.algebraic[36]*states[16]*states[17]
	    self.algebraic[63] = self.algebraic[59]/(self.algebraic[59]+self.algebraic[61])
	    self.algebraic[65] = 2.00000*(power(self.algebraic[64], 2.00000))*(power(self.algebraic[63], 2.00000))
	    self.algebraic[67] = (8.00000/3.00000)*(power(self.algebraic[64], 3.00000))*self.algebraic[63]
	    self.algebraic[69] = power(self.algebraic[64], 4.00000)
	    self.algebraic[70] = constants[93]*(states[0]-self.algebraic[56])*(self.algebraic[69]+self.algebraic[67]+self.algebraic[65])*states[18]
	    self.algebraic[71] = constants[94]*(states[0]-self.algebraic[56])*(0.600000*states[19]+0.400000*states[20])*states[21]
	    self.algebraic[72] = constants[32]*self.algebraic[39]*(power(states[22], 2.00000))*(0.900000*states[23]+0.100000)
	    self.algebraic[73] = constants[33]*self.algebraic[32]*(power(states[22], 2.00000))*(0.900000*states[23]+0.100000)
	    self.algebraic[74] = self.algebraic[73]+self.algebraic[72]
	    self.algebraic[75] = constants[35]*self.algebraic[39]*(power(states[24], 3.00000))*states[25]
	    self.algebraic[76] = constants[36]*self.algebraic[32]*(power(states[24], 3.00000))*states[25]
	    self.algebraic[77] = self.algebraic[76]+self.algebraic[75]
	    self.algebraic[101] = constants[9]*exp((-0.820000*constants[2]*states[0])/(constants[0]*constants[1]))
	    self.algebraic[102] = 1.00000/(1.00000+(power(constants[54]/self.algebraic[101], 1.06000))*(1.00000+power(constants[11]/constants[56], 1.12000)))
	    self.algebraic[98] = 1.00000/(1.00000+(power(constants[53]/states[1], 1.06000))*(1.00000+power(states[2]/constants[55], 1.12000)))
	    self.algebraic[100] = 0.370000/(1.00000+constants[57]/states[4])
	    self.algebraic[103] = constants[58]*constants[3]*1.00000*(self.algebraic[100]*self.algebraic[98]*states[27]-constants[50]*self.algebraic[102]*(1.00000-states[27]))
	    self.algebraic[94] = constants[48]*constants[3]*1.00000*(self.algebraic[92]*self.algebraic[90]*states[26]-self.algebraic[93]*constants[97]*(1.00000-states[26]))
	    self.algebraic[105] = self.algebraic[45]+self.algebraic[50]+self.algebraic[54]+self.algebraic[70]+self.algebraic[71]+self.algebraic[74]+self.algebraic[77]+self.algebraic[89]+self.algebraic[103]+self.algebraic[94]
	    self.rates[0] = -(self.algebraic[105]+self.algebraic[0])/constants[3]
	    self.algebraic[106] = self.algebraic[42]+self.algebraic[73]+self.algebraic[76]+self.algebraic[48]+self.algebraic[79]+self.algebraic[84]+3.00000*self.algebraic[103]+3.00000*self.algebraic[94]
	    self.rates[1] = -self.algebraic[106]/(constants[2]*constants[12])
	    self.algebraic[107] = (self.algebraic[70]+self.algebraic[71]+self.algebraic[75]+self.algebraic[87]+self.algebraic[72]+self.algebraic[44]+self.algebraic[49]+self.algebraic[78]+self.algebraic[83]+self.algebraic[81])-2.00000*self.algebraic[103]
	    self.rates[2] = -(self.algebraic[107]+self.algebraic[0])/(constants[2]*constants[12])
	    self.algebraic[104] = 1.00000/(1.00000+(power(constants[56]/constants[11], 1.12000))*(1.00000+power(self.algebraic[101]/constants[54], 1.06000)))
	    self.algebraic[109] = constants[50]*self.algebraic[102]+constants[52]*self.algebraic[104]
	    self.algebraic[99] = 1.00000/(1.00000+(power(constants[55]/states[2], 1.12000))*(1.00000+power(states[1]/constants[53], 1.06000)))
	    self.algebraic[111] = self.algebraic[100]*self.algebraic[98]+constants[51]*self.algebraic[99]
	    self.rates[27] = self.algebraic[109]*(1.00000-states[27])-self.algebraic[111]*states[27]
	    self.algebraic[110] = 1.00000/(1.00000+constants[63]/self.algebraic[29])
	    self.algebraic[108] = 1.00000/(1.00000+constants[62]/states[28])
	    self.algebraic[114] = 1.00000/(1.00000+constants[64]/states[4])
	    self.algebraic[115] = constants[65]*1.00000*(constants[59]*self.algebraic[108]*states[29]-self.algebraic[114]*self.algebraic[110]*(1.00000-states[29]))
	    self.algebraic[117] = -0.400000*states[34]*constants[75]
	    self.rates[4] = ((constants[15]*(constants[16]-states[4])+self.algebraic[117])-self.algebraic[103]/(constants[2]*constants[12]))+self.algebraic[115]/(4.00000*constants[2]*constants[12])
	    self.algebraic[113] = 1.00000-self.algebraic[110]
	    self.algebraic[116] = self.algebraic[114]*self.algebraic[110]+constants[61]*self.algebraic[113]
	    self.algebraic[112] = 1.00000-self.algebraic[108]
	    self.algebraic[118] = constants[59]*self.algebraic[108]+constants[60]*self.algebraic[112]
	    self.rates[29] = self.algebraic[116]*(1.00000-states[29])-self.algebraic[118]*states[29]
	    self.algebraic[52] = self.algebraic[51]*self.algebraic[46]
	    self.algebraic[119] = 280000.*(power(self.algebraic[29]/1.00000, 2.00000))+constants[68]*self.algebraic[52]
	    self.algebraic[121] = (constants[73]-states[32])+constants[74]
	    self.algebraic[122] = constants[74]*states[32]
	    self.algebraic[123] = (power(power(self.algebraic[121], 2.00000)+4.00000*self.algebraic[122], 1.0/2)-self.algebraic[121])/2.00000
	    self.algebraic[124] = 0.0800000/(1.00000+0.360000/self.algebraic[123])
	    self.rates[30] = states[31]*self.algebraic[119]-states[30]*self.algebraic[124]
	    self.algebraic[125] = 0.000377000*(power(self.algebraic[123]/1.00000, 2.00000))
	    self.rates[31] = self.algebraic[125]*(1.00000-(states[30]+states[31]))-(self.algebraic[119]+constants[67])*states[31]
	    self.algebraic[126] = constants[66]*(self.algebraic[123]-self.algebraic[29])*states[30]
	    self.algebraic[127] = constants[69]*(states[28]-self.algebraic[123])
	    self.rates[32] = (self.algebraic[127]-self.algebraic[126])/(2.00000*constants[2]*constants[71])
	    self.algebraic[120] = constants[70]*(states[28]-self.algebraic[29])
	    self.rates[28] = ((-self.algebraic[115]-self.algebraic[127])-self.algebraic[120])/(2.00000*constants[2]*constants[72])
	    self.algebraic[134] = constants[77]*states[33]*constants[98]-constants[82]*states[34]
	    self.algebraic[130] = ((1.00000-states[33])-states[34])-states[35]
	    self.algebraic[132] = constants[76]*self.algebraic[29]*self.algebraic[130]-constants[81]*states[33]
	    self.rates[33] = self.algebraic[132]-self.algebraic[134]
	    self.algebraic[128] = constants[85]-states[36]
	    self.rates[36] = constants[90]*(self.algebraic[128]-constants[91])
	    self.algebraic[136] = constants[78]*states[34]-constants[83]*states[35]*self.algebraic[29]
	    self.algebraic[139] = constants[80]*(power(self.rates[36], 2.00000))*states[34]
	    self.rates[34] = (self.algebraic[134]-self.algebraic[136])-self.algebraic[139]
	    self.algebraic[137] = constants[79]*states[35]
	    self.algebraic[138] = constants[80]*(power(self.rates[36], 2.00000))*states[35]
	    self.rates[35] = (self.algebraic[136]-self.algebraic[137])-self.algebraic[138]
	    self.algebraic[96] = (self.algebraic[47]+self.algebraic[54]+self.algebraic[88])-2.00000*self.algebraic[94]
	    self.algebraic[140] = constants[75]*((self.algebraic[139]+self.algebraic[136])-self.algebraic[132])
	    self.rates[3] = -(((self.algebraic[96]-self.algebraic[115])-self.algebraic[126])-self.algebraic[120])/(2.00000*constants[2]*constants[12])+self.algebraic[140]
	    # print "self.rates is : ",self.rates
	    return(self.rates)

	def computeAlgebraic(self,voi,constants,states):
	    self.algebraic[2] = 1.00000/(9.00000e+09*exp(states[0]/5.00000)+8000.00*exp(states[0]/100.000))
	    self.algebraic[15] = 1.00000/(0.0140000*exp(-states[0]/5.00000)+4000.00*exp(-states[0]/100.000))
	    self.algebraic[4] = 1.00000/(250000.*exp(states[0]/9.00000)+58.0000*exp(states[0]/65.0000))
	    self.algebraic[17] = 1.00000/(1800.00*exp(-states[0]/14.0000)+66.0000*exp(-states[0]/65.0000))
	    self.algebraic[5] = 1.00000/(0.0190000*exp(-states[0]/5.60000)+0.820000*exp(-states[0]/250.000))
	    self.algebraic[18] = 1.00000/(40.0000*exp(states[0]/6.30000)+1.50000*exp(states[0]/10000.0))
	    self.algebraic[6] = 1.00000/(62000.0*exp(states[0]/10.1000)+30.0000*exp(states[0]/3000.00))
	    self.algebraic[19] = 1.00000/(0.000600000*exp(-states[0]/6.70000)+1.20000*exp(-states[0]/25.0000))
	    self.algebraic[7] = 1.00000/(20.0000*exp(-states[0]/11.5000)+5.00000*exp(-states[0]/300.000))
	    self.algebraic[20] = 1.00000/(160.000*exp(states[0]/28.0000)+200.000*exp(states[0]/1000.00))+1.00000/(2500.00*exp(states[0]/20.0000))
	    self.algebraic[8] = 1.00000/(200.000*exp(-states[0]/13.0000)+20.0000*exp(-states[0]/300.000))
	    self.algebraic[21] = 1.00000/(1600.00*exp(states[0]/28.0000)+2000.00*exp(states[0]/1000.00))+1.00000/(10000.0*exp(states[0]/20.0000))
	    self.algebraic[9] = 1.00000/(10.0000*exp(states[0]/17.0000)+2.50000*exp(states[0]/300.000))
	    self.algebraic[22] = 1.00000/(0.350000*exp(-states[0]/17.0000)+2.00000*exp(-states[0]/150.000))
	    self.algebraic[10] = 1.00000/(85.0000*exp(-states[0]/10.5000)+370.000*exp(-states[0]/62.0000))
	    self.algebraic[23] = 1.00000/(1450.00*exp(states[0]/20.0000)+260.000*exp(states[0]/100.000))
	    self.algebraic[11] = 1.00000/(11.0000*exp(-states[0]/28.0000)+0.200000*exp(-states[0]/400.000))
	    self.algebraic[24] = 1.00000/(4.40000*exp(states[0]/16.0000)+0.200000*exp(states[0]/500.000))
	    self.algebraic[12] = (0.00380000*exp(-(states[0]+13.5000)/11.3000))/(1.00000+0.0513350*exp(-(states[0]+13.5000)/11.3000))
	    self.algebraic[25] = (0.00380000*exp((states[0]+13.5000)/11.3000))/(1.00000+0.0670830*exp((states[0]+13.5000)/11.3000))
	    self.algebraic[16] = 1.00000/(0.270000*exp(-states[0]/5.90000)+1.50000*exp(-states[0]/65.0000))
	    self.algebraic[28] = 1.00000/(480.000*exp(states[0]/7.00000)+2.20000*exp(states[0]/65.0000))
	    self.algebraic[14] = 1.00000/(0.102700*exp(-states[0]/8.00000)+0.250000*exp(-states[0]/50.0000))
	    self.algebraic[27] = 1.00000/(26.0000*exp(states[0]/17.0000)+0.0200000*exp(states[0]/800.000))
	    self.algebraic[30] = 1.00000/(0.800000*exp(-states[0]/400.000))
	    self.algebraic[3] = ((1.00000-states[9])-states[13])-states[14]
	    self.algebraic[31] = 1.00000/(0.00180000*exp(-states[0]/7.40000)+2.00000*exp(-states[0]/100.000))
	    self.algebraic[34] = 1.00000/(2.20000e+06*exp(states[0]/7.40000)+11.0000*exp(states[0]/100.000))
	    self.algebraic[13] = (constants[13]-states[3])+constants[14]
	    self.algebraic[26] = constants[14]*states[3]
	    self.algebraic[29] = (power(power(self.algebraic[13], 2.00000)+4.00000*self.algebraic[26], 1.0/2)-self.algebraic[13])/2.00000
	    self.algebraic[35] = 3.70000*self.algebraic[29]
	    self.algebraic[1] = ((1.00000-states[7])-states[5])-states[8]
	    self.algebraic[37] = 1.00000/(0.000102700*exp(-states[0]/8.00000)+5.00000*exp(-states[0]/400.000))
	    self.algebraic[33] = 1.00000/(1300.00*exp(states[0]/20.0000)+0.0400000*exp(states[0]/800.000))
	    self.algebraic[38] = 0.0400000/(1.00000+(constants[21]*self.algebraic[28]*self.algebraic[31])/(constants[20]*self.algebraic[16]*self.algebraic[34]))
	    self.algebraic[41] = 0.0400000-self.algebraic[38]
	    self.algebraic[40] = 0.0100000/(1.00000+(constants[18]*self.algebraic[27]*self.algebraic[37])/(self.algebraic[30]*self.algebraic[14]*self.algebraic[33]))
	    self.algebraic[43] = 0.0100000-self.algebraic[40]
	    self.algebraic[36] = self.custom_piecewise([equal(states[0] , 0.00000), -constants[10] , True, (((2.00000*constants[2]*states[0])/(constants[0]*constants[1]))*(self.algebraic[29]-constants[10]*exp((-2.00000*constants[2]*states[0])/(constants[0]*constants[1]))))/(1.00000-exp((-2.00000*constants[2]*states[0])/(constants[0]*constants[1])))])
	    self.algebraic[51] = 0.0676000*self.algebraic[36]
	    self.algebraic[53] = self.algebraic[29]-0.300000*self.algebraic[51]
	    self.algebraic[55] = self.algebraic[53]*states[9]
	    self.algebraic[57] = self.algebraic[55]+self.algebraic[29]*(1.00000-states[9])
	    self.algebraic[60] = constants[26]*self.algebraic[57]
	    self.algebraic[62] = ((1.00000-states[15])-states[10])-states[11]
	    self.algebraic[56] = ((constants[0]*constants[1])/constants[2])*log(constants[11]/states[2])
	    self.algebraic[66] = 1.00000/(8000.00*exp(((states[0]-self.algebraic[56])-97.0000)/8.50000)+7.00000*exp(((states[0]-self.algebraic[56])-97.0000)/300.000))
	    self.algebraic[59] = (0.750000*exp(0.0350000*((states[0]-self.algebraic[56])-10.0000)))/(1.00000+exp(0.0150000*((states[0]-self.algebraic[56])-140.000)))
	    self.algebraic[61] = (3.00000*exp(-0.0480000*((states[0]-self.algebraic[56])-10.0000))*(1.00000+exp(0.0640000*((states[0]-self.algebraic[56])-38.0000))))/(1.00000+exp(0.0300000*((states[0]-self.algebraic[56])-70.0000)))
	    self.algebraic[64] = self.algebraic[61]/(self.algebraic[59]+self.algebraic[61])
	    self.algebraic[68] = ((power(self.algebraic[64], 4.00000))*1.00000)/(0.000140000*exp(-((states[0]-self.algebraic[56])-97.0000)/9.10000)+0.200000*exp(-((states[0]-self.algebraic[56])-97.0000)/500.000))
	    self.algebraic[93] = 1.00000*exp(((constants[49]-1.00000)*constants[2]*states[0])/(constants[0]*constants[1]))
	    self.algebraic[95] = self.algebraic[93]*constants[97]+constants[43]*constants[100]
	    self.algebraic[90] = 1.00000/(1.00000+(power(constants[44]/states[1], 3.00000))*(1.00000+self.algebraic[29]/constants[46]))
	    self.algebraic[91] = 1.00000/(1.00000+(constants[46]/self.algebraic[29])*(1.00000+power(states[1]/constants[44], 3.00000)))
	    self.algebraic[92] = 1.00000*exp((constants[49]*constants[2]*states[0])/(constants[0]*constants[1]))
	    self.algebraic[97] = self.algebraic[92]*self.algebraic[90]+constants[42]*self.algebraic[91]
	    self.algebraic[0] = self.custom_piecewise([greater_equal(voi , constants[4]) & less_equal(voi , constants[5]) & less_equal((voi-constants[4])-floor((voi-constants[4])/constants[6])*constants[6] , constants[7]), constants[8] , True, 0.00000])
	    self.algebraic[32] = self.custom_piecewise([equal(states[0] , 0.00000), -constants[9] , True, (((constants[2]*states[0])/(constants[0]*constants[1]))*(states[1]-constants[9]*exp((-constants[2]*states[0])/(constants[0]*constants[1]))))/(1.00000-exp((-constants[2]*states[0])/(constants[0]*constants[1])))])
	    self.algebraic[79] = constants[37]*self.algebraic[32]
	    self.algebraic[39] = self.custom_piecewise([equal(states[0] , 0.00000), states[2] , True, (((constants[2]*states[0])/(constants[0]*constants[1]))*(states[2]-constants[11]*exp((-constants[2]*states[0])/(constants[0]*constants[1]))))/(1.00000-exp((-constants[2]*states[0])/(constants[0]*constants[1])))])
	    self.algebraic[78] = 0.400000*constants[37]*self.algebraic[39]
	    self.algebraic[80] = self.algebraic[78]+self.algebraic[79]
	    self.algebraic[88] = constants[41]*self.algebraic[36]
	    self.algebraic[81] = self.custom_piecewise([equal(states[0] , -3.00000), constants[95]*self.algebraic[39]*13.0077 , True, (constants[95]*self.algebraic[39]*(states[0]+3.00000))/(1.00000-exp(-(states[0]+3.00000)/13.0000))])
	    self.algebraic[82] = 1.00000/(1.00000+power(0.00120000/self.algebraic[29], 3.00000))
	    self.algebraic[84] = constants[38]*self.algebraic[32]*self.algebraic[82]
	    self.algebraic[83] = constants[38]*self.algebraic[39]*self.algebraic[82]
	    self.algebraic[85] = self.algebraic[83]+self.algebraic[84]
	    self.algebraic[86] = 0.800000/(1.00000+power(states[4]/0.100000, 2.00000))
	    self.algebraic[87] = constants[96]*(states[0]-self.algebraic[56])*self.algebraic[86]
	    self.algebraic[89] = self.algebraic[80]+self.algebraic[88]+self.algebraic[81]+self.algebraic[85]+self.algebraic[87]
	    self.algebraic[42] = constants[17]*self.algebraic[32]*states[5]*states[6]
	    self.algebraic[44] = 0.100000*constants[17]*self.algebraic[39]*states[5]*states[6]
	    self.algebraic[45] = self.algebraic[42]+self.algebraic[44]
	    self.algebraic[46] = (states[9]*(states[10]+states[11])*states[12])/(1.00000+power(1.40000/states[4], 3.00000))
	    self.algebraic[48] = 1.85000e-05*constants[19]*self.algebraic[32]*self.algebraic[46]
	    self.algebraic[49] = 0.000365000*constants[19]*self.algebraic[39]*self.algebraic[46]
	    self.algebraic[47] = constants[19]*self.algebraic[36]*self.algebraic[46]
	    self.algebraic[50] = self.algebraic[48]+self.algebraic[47]+self.algebraic[49]
	    self.algebraic[54] = constants[29]*self.algebraic[36]*states[16]*states[17]
	    self.algebraic[63] = self.algebraic[59]/(self.algebraic[59]+self.algebraic[61])
	    self.algebraic[65] = 2.00000*(power(self.algebraic[64], 2.00000))*(power(self.algebraic[63], 2.00000))
	    self.algebraic[67] = (8.00000/3.00000)*(power(self.algebraic[64], 3.00000))*self.algebraic[63]
	    self.algebraic[69] = power(self.algebraic[64], 4.00000)
	    self.algebraic[70] = constants[93]*(states[0]-self.algebraic[56])*(self.algebraic[69]+self.algebraic[67]+self.algebraic[65])*states[18]
	    self.algebraic[71] = constants[94]*(states[0]-self.algebraic[56])*(0.600000*states[19]+0.400000*states[20])*states[21]
	    self.algebraic[72] = constants[32]*self.algebraic[39]*(power(states[22], 2.00000))*(0.900000*states[23]+0.100000)
	    self.algebraic[73] = constants[33]*self.algebraic[32]*(power(states[22], 2.00000))*(0.900000*states[23]+0.100000)
	    self.algebraic[74] = self.algebraic[73]+self.algebraic[72]
	    self.algebraic[75] = constants[35]*self.algebraic[39]*(power(states[24], 3.00000))*states[25]
	    self.algebraic[76] = constants[36]*self.algebraic[32]*(power(states[24], 3.00000))*states[25]
	    self.algebraic[77] = self.algebraic[76]+self.algebraic[75]
	    self.algebraic[101] = constants[9]*exp((-0.820000*constants[2]*states[0])/(constants[0]*constants[1]))
	    self.algebraic[102] = 1.00000/(1.00000+(power(constants[54]/self.algebraic[101], 1.06000))*(1.00000+power(constants[11]/constants[56], 1.12000)))
	    self.algebraic[98] = 1.00000/(1.00000+(power(constants[53]/states[1], 1.06000))*(1.00000+power(states[2]/constants[55], 1.12000)))
	    self.algebraic[100] = 0.370000/(1.00000+constants[57]/states[4])
	    self.algebraic[103] = constants[58]*constants[3]*1.00000*(self.algebraic[100]*self.algebraic[98]*states[27]-constants[50]*self.algebraic[102]*(1.00000-states[27]))
	    self.algebraic[94] = constants[48]*constants[3]*1.00000*(self.algebraic[92]*self.algebraic[90]*states[26]-self.algebraic[93]*constants[97]*(1.00000-states[26]))
	    self.algebraic[105] = self.algebraic[45]+self.algebraic[50]+self.algebraic[54]+self.algebraic[70]+self.algebraic[71]+self.algebraic[74]+self.algebraic[77]+self.algebraic[89]+self.algebraic[103]+self.algebraic[94]
	    self.algebraic[106] = self.algebraic[42]+self.algebraic[73]+self.algebraic[76]+self.algebraic[48]+self.algebraic[79]+self.algebraic[84]+3.00000*self.algebraic[103]+3.00000*self.algebraic[94]
	    self.algebraic[107] = (self.algebraic[70]+self.algebraic[71]+self.algebraic[75]+self.algebraic[87]+self.algebraic[72]+self.algebraic[44]+self.algebraic[49]+self.algebraic[78]+self.algebraic[83]+self.algebraic[81])-2.00000*self.algebraic[103]
	    self.algebraic[104] = 1.00000/(1.00000+(power(constants[56]/constants[11], 1.12000))*(1.00000+power(self.algebraic[101]/constants[54], 1.06000)))
	    self.algebraic[109] = constants[50]*self.algebraic[102]+constants[52]*self.algebraic[104]
	    self.algebraic[99] = 1.00000/(1.00000+(power(constants[55]/states[2], 1.12000))*(1.00000+power(states[1]/constants[53], 1.06000)))
	    self.algebraic[111] = self.algebraic[100]*self.algebraic[98]+constants[51]*self.algebraic[99]
	    self.algebraic[110] = 1.00000/(1.00000+constants[63]/self.algebraic[29])
	    self.algebraic[108] = 1.00000/(1.00000+constants[62]/states[28])
	    self.algebraic[114] = 1.00000/(1.00000+constants[64]/states[4])
	    self.algebraic[115] = constants[65]*1.00000*(constants[59]*self.algebraic[108]*states[29]-self.algebraic[114]*self.algebraic[110]*(1.00000-states[29]))
	    self.algebraic[117] = -0.400000*states[34]*constants[75]
	    self.algebraic[113] = 1.00000-self.algebraic[110]
	    self.algebraic[116] = self.algebraic[114]*self.algebraic[110]+constants[61]*self.algebraic[113]
	    self.algebraic[112] = 1.00000-self.algebraic[108]
	    self.algebraic[118] = constants[59]*self.algebraic[108]+constants[60]*self.algebraic[112]
	    self.algebraic[52] = self.algebraic[51]*self.algebraic[46]
	    self.algebraic[119] = 280000.*(power(self.algebraic[29]/1.00000, 2.00000))+constants[68]*self.algebraic[52]
	    self.algebraic[121] = (constants[73]-states[32])+constants[74]
	    self.algebraic[122] = constants[74]*states[32]
	    self.algebraic[123] = (power(power(self.algebraic[121], 2.00000)+4.00000*self.algebraic[122], 1.0/2)-self.algebraic[121])/2.00000
	    self.algebraic[124] = 0.0800000/(1.00000+0.360000/self.algebraic[123])
	    self.algebraic[125] = 0.000377000*(power(self.algebraic[123]/1.00000, 2.00000))
	    self.algebraic[126] = constants[66]*(self.algebraic[123]-self.algebraic[29])*states[30]
	    self.algebraic[127] = constants[69]*(states[28]-self.algebraic[123])
	    self.algebraic[120] = constants[70]*(states[28]-self.algebraic[29])
	    self.algebraic[134] = constants[77]*states[33]*constants[98]-constants[82]*states[34]
	    self.algebraic[130] = ((1.00000-states[33])-states[34])-states[35]
	    self.algebraic[132] = constants[76]*self.algebraic[29]*self.algebraic[130]-constants[81]*states[33]
	    self.algebraic[128] = constants[85]-states[36]
	    self.algebraic[136] = constants[78]*states[34]-constants[83]*states[35]*self.algebraic[29]
	    self.algebraic[139] = constants[80]*(power(self.rates[36], 2.00000))*states[34]
	    self.algebraic[137] = constants[79]*states[35]
	    self.algebraic[138] = constants[80]*(power(self.rates[36], 2.00000))*states[35]
	    self.algebraic[96] = (self.algebraic[47]+self.algebraic[54]+self.algebraic[88])-2.00000*self.algebraic[94]
	    self.algebraic[140] = constants[75]*((self.algebraic[139]+self.algebraic[136])-self.algebraic[132])
	    self.algebraic[58] = constants[28]*self.algebraic[55]
	    self.algebraic[129] = constants[75]*(states[34]+states[35])
	    self.algebraic[131] = constants[89]*self.algebraic[129]
	    self.algebraic[133] = self.algebraic[131]*self.algebraic[128]
	    self.algebraic[135] = -constants[99]+self.algebraic[133]
	    # return self.algebraic

	def custom_piecewise(self, cases):
	    """Compute result of a piecewise function"""
	    return select(cases[0::2],cases[1::2])

	def step_model(self, constants, states, circumference):
		"""Solve model with ODE solver"""

	    # # Construct ODE object to solve
	    # self.r = ode(computeRates)
	    # self.r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
	    # self.r.set_initial_value(init_self.states, voi[0])
	    # self.r.set_f_params(constants)

	    # # Solve model
	    # states = array([[0.0] * len(voi)] * self.sizeself.states)
	    # states[:,0] = init_self.states
	    # for (i,t) in enumerate(voi[1:]):
		# states = [-1.0,-1.0]
		if self.r.successful():
			savedStates = self.r.y
			savedTime = self.r.t
			for i in range(0,3):
				self.r.set_f_params(constants)
				self.r.integrate(self.time)
				self.halfSarcomereLength =  circumference / self.numberOfMuscleUnitsAroundCircumference * 1000 # *1000 is to convert to um # originally this was: constants[85] = 0.9623799975411884
				# self.halfSarcomereLength = halfSarcomereLength = 0.9623799975411884
				(init_states, constants) = self.initConsts(self.halfSarcomereLength)
				self.r.set_initial_value(savedStates, savedTime)

			# print "self states 1:", states[1]
			# print "doing self.time: ", self.time
			self.r.set_f_params(constants)
			self.r.integrate(self.time)
			states = self.r.y
        # else:
        #     break

	    # Compute self.algebraic variables
		# self.computeAlgebraic(self.time,constants)
	    # return (voi, states, self.algebraic)
		return states

	def __init__(self):
		self.stepIndex = 0
		self.tensionHistory = array([0.0]*200000)
		self.sarcomereHalfLengthHistory = array([0.0]*200000)
		self.tensionHistory = array([0.0]*200000)
		self.imposedPressureHistory = array([0.0]*200000)

		self.customOutputArray1 = array([0.0]*200000)
		self.customOutputArray2 = array([0.0]*200000)
		# self.m_periodicTime = 0.0; #\todo think about this for restarts!
		# self.m_timeToMaximumElastance = 0.2782;
		# self.m_timeToRelax = 0.1391;
		self.m_minimumElastance = 4.10246e-3;
		self.m_maximumElastance = 3.0827e-1/35.0;
		
		# self.m_heartPeriod = 0.86;
		# Size of variable arrays:
		self.sizeAlgebraic = 141
		self.sizeStates = 37
		self.sizeConstants = 101
		self.rates = []

		self.time = 0.0
		self.stimPeriod = 800 #ms
		self.stimStartTime = 50
		self.stimDuration = 2

		self.halfSarcomereLength = 0.9623799975411884

		self.algebraic = array([0.0] * self.sizeAlgebraic)

		# Initialise constants and state variables
		originalCellMLValueOfConstants85 = 0.9623799975411884
		(init_states, constants) = self.initConsts(originalCellMLValueOfConstants85)

		loadInitialStatesFromFile = True
		if loadInitialStatesFromFile:
			self.states = loadtxt('eightHundredMsInitialConditions_sm03.dat')


	    # Construct ODE object to solve
		self.r = ode(self.computeRates)
		self.r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
		self.r.set_initial_value(self.states, self.time)
		self.r.set_f_params(constants)


		# for i in range(0,5000):
		# 	elastance = self.updateControl(-1, 0.2, -1, -1)
		# 	print elastance







		# for i in range(0,1000):
		# 	self.time = self.time + 0.2
		# 	self.states = self.step_model(constants,self.states)

		# 	self.computeAlgebraic(self.time, constants, self.states)
		# 	elastance = self.algebraic[135] * self.m_maximumElastance
		# 	print elastance


	# def solve_model():
	#     """Solve model with ODE solver"""
	    

	#     # Set timespan to solve over
	#     simulationLengthInMilliseconds = 2000
	#     dtInMilliseconds = 0.2
	#     totalSteps = simulationLengthInMilliseconds / dtInMilliseconds
	#     voi = linspace(0, simulationLengthInMilliseconds, totalSteps)

	#     # Construct ODE object to solve
	#     self.r = ode(computeRates)
	#     self.r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
	#     self.r.set_initial_value(init_self.states, voi[0])
	#     self.r.set_f_params(constants)

	#     # Solve model
	#     states = array([[0.0] * len(voi)] * self.sizeself.states)
	#     states[:,0] = init_self.states
	#     for (i,t) in enumerate(voi[1:]):
	#         if self.r.successful():
	#             self.r.integrate(t)
	#             states[:,i+1] = self.r.y
	#         else:
	#             break

	#     # Compute self.algebraic variables
	#     self.algebraic = computeAlgebraic(constants, states, voi)
	#     return (voi, states, self.algebraic)

	# def plot_model(voi, states, self.algebraic):
	#     """Plot variables against variable of integration"""
	#     import pylab
	#     (legend_self.states, legend_algebraic, legend_voi, legend_constants) = createLegends()
	#     pylab.figure(1)
	#     saveself.states = False
	#     if saveself.states:
	#         savetxt('eightHundredMsInitial_sm03',states)
	#     # pylab.plot(voi,vstack((states,self.algebraic)).T)
	#     pylab.plot(voi,vstack((self.algebraic[133:135])).T)
	#     pylab.xlabel(legend_voi)
	#     # pylab.legend(legend_self.states + legend_algebraic, loc='best')
	#     pylab.legend(legend_algebraic[133:135], loc='best')
	#     pylab.show()

	# if __name__ == "__main__":
	#     (voi, states, self.algebraic) = solve_model()
	    # plot_model(voi, states, self.algebraic)

# if __name__ == "__main__":
# 	a= thinShellPressureController()
# 	a.updateControl(1, 0.01, [], [], [1,2,3,4,5,6,7,8,9])
