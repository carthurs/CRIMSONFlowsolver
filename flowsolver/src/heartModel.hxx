#ifndef HEARTMODEL_HXX_
#define HEARTMODEL_HXX_

class heart
{
public:
	heart(int surfaceIndex_in)
	{
		
	}
private:
	double  patrial;                     // atrial pressure
	double  edv;                         //
	double  rmv;                         // mitral valve resistance
	double  lmv;						 // mitral valve inductance
	double  qmv_n;                       // mitral valve flow at t=n and t=n+1         
	double  rav;                         // aortic valve resistance
	double  lav;						 // aoritc valve inductance
	double  vlv_n;                       // left ventricular volume at t=n and t=n+1
	double  plv_n;                       // left ventricular pressure at t=n and t=n+1         
	double  vulv;                        // left ventricular unstressed volume
	double  pla_n;                       // left atrial pressure at t=n and t=n+1         
	double  vula;                        // left atrial unstressed volume
	double  emax;	                     // maximum myocardial elastance
	double  emin;						 // minimum myocardial elastance
	double  period;                      // heart period
	double  tmax;                        // 
	double  trelax;                      // 
	double  kelv;                        // defined value for elastance varying resistance 
	std::pair<double,double>  vlv_coeff; // vlv coefficients
	double  activationtime;              // elastance activation time
	std::vector<double>  vlv_hist; 	   			     // vlv history
	std::vector<double>  plv_hist; 				     // plv history         
	std::vector<double>  qmv_hist; 				     // qmv history       
	std::vector<double>  elv_hist; 			    	 // elv history            
	std::vector<std::vector<double>>  eval_hist;  	 // eval history            
	std::vector<double>  stab_hist;  			     // stabilisation pressure history
	bool inputHRandSP;                // Flag to indicate use of HeartRate.dat and SystolicPressure.dat input files
	type(linkedlist), pointer :: heartRateList // prescribed heart rate data, when we want to input this from data
	type(linkedlist), pointer :: systolicPressureList // prescribed peak systolic pressure array, for when we want to VERY APPROXIMATELY prescribe this from input data
	double totalTimeOfAllCompletedPeriods; // provides the time at which the present beat started
	bool avopen = false;             // aortic valve logical              
	type(timedata) :: elv_senzaki
	type(timedata) :: elv_input
	bool input_elastance;
	double evalval;
	int evalcount;
	*double sPress = 0;
	bool backflow = false;                 // backflow logical
	double  m_backflow;                    // backflow magnitude
	double  s_backflow;                    // backflow steepness 
	double  c_backflow;                    // backflow closure 
	double  t_backflow;                    // time spent in backflow         
	double  max_backflow = 39.5e-3;        // maximum backflow time = 39.5 ms from Leyh et al. Circulation, 1999, 100, 2153-2160          
	std::vector<double> act_hist;				       // activation history
	bool ibackflow = false;             //
	int timestepsSinceLastHeartPeriodReset; // counter to store how far through the current heart beat we are
	int newBeatJustStarted;            // This should only be set to 1 on the first timestep of the current beat, starting from the /second/ beat of the simulation. Otherwise, zero.
	character(len=50) :: ahistfilename     // activation time, normalised by period [0-1]
	character(len=50) :: avarsformat       // 

	void initialise_hrt();
	void initxvars_hrt();
	void updxvars_hrt();
	void setimplicitcoeff_hrt();
	void iterate_hrt();
	void writexvars_hrt();
	void isavopen();
	void getelastance_senzaki();
	void setpresspntrs_hrt();
	void update_activationtime_hrt();
	void write_activation_history_hrt();
	void read_activation_history_hrt();
	void set_sPress();
	void assign_ptrs_ext_hrt();
};

#endif