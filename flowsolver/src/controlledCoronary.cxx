#include "controlledCoronary.hxx"
#include <boost/math/special_functions/round.hpp>

// Statics
int controlledCoronary::numberOfInitialisedCoronaries = 0;

std::pair<double,double> controlledCoronary::computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt)
{
    // Get the intramyocardial pressure for the previous time-step and previous-previous time-step
    // note that these should both really be evaluated one time-step later, but to avoid the need
    // for a simultaneous solve of the heart and coronaries, we do it this way for now
    // \todo fix this
    if ((timestepNumber == int(0)) && (timestepNumber == int(1))) // treat case with no known IM pressure yet
    {
        P_IM_mid_lasttimestep = 0.0;//5000.0; // \todo find a better way of doing this; maybe input this value from file...
        P_IM_mid = 0.0;//5000.0; // ... or set it based on the aortic valve state at simulation start
    }
    else if (timestepNumber == int(2)) // treat case where only one IM pressure history point is known
    {
        P_IM_mid_lasttimestep = 0.0;
        P_IM_mid = 0.0;
        std::cout << "heart model needed for this!" << std::endl;
        // P_IM_mid_lasttimestep = this%intramyocardialPressureToLVScaling * hrt%plv_hist(timestepNumber-1)
        // P_IM_mid = this%intramyocardialPressureToLVScaling * hrt%plv_hist(timestepNumber)
    }
    else // get the previous intramyocardial pressure in the case where we have enough doata for this (see comment before start of "if" block)
    {
        P_IM_mid_lasttimestep = 0.0;
        P_IM_mid = 0.0;
        std::cout << "heart model needed for this!" << std::endl;
        // P_IM_mid_lasttimestep = this%intramyocardialPressureToLVScaling * hrt%plv_hist(timestepNumber-1) //\todo check these actually exist on first iteration
        // // Get IM pressure for now (this will be adjusted in a moment if we're on, according to alpha in gen alpha method)
        // P_IM_mid = this%intramyocardialPressureToLVScaling * hrt%plv_hist(timestepNumber) //\todo check these actually exist on first iteration
    }

    // we need these a couple of times to compute coeff, so precompute them for convenience.
    double m11 = 1.0 + resistanceNearAorta*complianceNearAorta/alfi_delt;
    double m12 = resistanceNearAorta * (1.0/distalResistance + intramyocardialCompliance/alfi_delt);
    double m22 = midResistance * (intramyocardialCompliance/alfi_delt + 1.0/distalResistance) + 1.0;
    double temporary_variable = m22/(m11*m22-m22+m12);

    // The actual differential equation "solve" - here we obtain the operators
    // which describe the boundary condtion from the ODE.
    double temp1 = (m11 + m12/m22) * resistanceNearAorta * temporary_variable;
    double temp2 = ((capacitorNearAortaTopPressure*complianceNearAorta + intramyocardialCapacitorTopPressure*intramyocardialCompliance + P_IM_mid*intramyocardialCompliance - P_IM_mid_lasttimestep
                   *intramyocardialCompliance)*resistanceNearAorta/alfi_delt - m12/m22*midResistance*intramyocardialCompliance/alfi_delt
                   *(intramyocardialCapacitorTopPressure+P_IM_mid-P_IM_mid_lasttimestep))*temporary_variable;


    std::pair<double,double> returnCoeffs;
    returnCoeffs.first = temp1;
    returnCoeffs.second = temp2;

    std::cout << "Coefficient 1: " << temp1 << std::endl;
    std::cout << "Coefficient 2: " << temp2 << std::endl;

    return returnCoeffs;
}

void controlledCoronary::initialiseModel()
{
	// allocate(this%deltaMVO2(this%numberOfControlledCoronaryModelSurfaces))
	// allocate(this%deltaMVO2PerTimestepOverPreviousBeat(this%numberOfControlledCoronaryModelSurfaces)) // There's some redundancy between this and deltaMVO2 at the time of writing, but the separation is likely to be useful in future modifications
	// allocate(this%P_a(this%numberOfControlledCoronaryModelSurfaces))

	MVO2computedAtEndOfBeatPreviousToPreviousBeat = 0.0;
	MVO2computedAtEndOfPreviousBeat = 0.0;


	deltaMVO2PerTimestepOverPreviousBeat = 0.0;

    // these values should never be used as they are re-assigned before the first time they're read, but set zero anyway
	P_IM_mid = 0.0;
	P_IM_mid_lasttimestep = 0.0;

    currentMyocardialHungerSignal = 0.0;
    hungerDelta = 0.0;

	// These are currently set to match the RCR test. Probably it'd be better
    // to set these manually!
	// capacitorNearAortaTopPressure_acceptedAtEndOfLastTimestep = 10585.0000000000;//10000.0; //\todo make these settable for mmgs and cmgs
    capacitorNearAortaTopPressure = capacitorNearAortaTopPressure_acceptedAtEndOfLastTimestep;

	// intramyocardialCapacitorTopPressure_acceptedAtEndOfLastTimestep = 6.080839425964667e+03;//10000.0;
    intramyocardialCapacitorTopPressure = intramyocardialCapacitorTopPressure_acceptedAtEndOfLastTimestep;


	if (thisIsARestartedSimulation)
    {
		// Initialise the pressure using the value from the PHistRCR.dat.
		std::cout << "this is all wrong... restart not set up yet!" << std::endl;
		std::exit(1);
		// pressure_n = (boundaryConditionManager::Instance()->PHistReader)->getReadFileData(indexOfThisCoronary+1,timdat.lstep);
    }
    else
    {
		pressure_n = *pressure_n_ptr;

        std::cout << "press pointer " << pressure_n_ptr << " " << *pressure_n_ptr << std::endl;
        std::cout << "flow ptr " << flow_n_ptr << " " << *flow_n_ptr << std::endl;
    }
	
	intramyocardialPressureToLVScaling = 1.0; // \todo try other values here (0.4 in MATLAB)

    // iround is a to-integer rounding function
	O2supplyDemandHistoryWindowLength_timesteps = boost::math::iround<double>(O2supplyDemandHistoryWindowLength_seconds/delt);

	// Ensure that push_back on this vector does not cause big in-memory copies by reserving (empty)
	// memory for the vector to grow into.
	O2supplyDemandHistory.reserve(inpdat.nstep[0] + O2supplyDemandHistoryWindowLength_timesteps + 1);
	// set the history of the supply-demand discrepancy at each outlet (note that this array is long enough
	//  for the entire simulation + a fake "negative time" history of zeros; we set this fake history now:
	//\todo worry about restart issues here.
	O2supplyDemandHistory.resize(O2supplyDemandHistoryWindowLength_timesteps,0.0);
	
	// initialise the moving average O2 supply/demand history to zero
	currentMyocardialHungerSignal = 0.0; // \todo have a way of setting the initial discrepancy history to be non-zero
	
	// set oneOverTotalResistance to be one over the initial total effective resistance of the LPN
	oneOverTotalResistance = 1.0/(resistanceNearAorta + midResistance + distalResistance);

	arterialO2VolumeProportion = 7.0/40.0; // percent \todo this could be dynamic!

// ----------------------NEW CORONARY RESTART CODE COMPONENT----------------------
         // // Test whether this is a restart; if so, overwrite the necessary values with those from disk:
         // inquire(file="coronary_restart.dat", exist=file_exists)
         // if (file_exists) then
         // // \todo confirm this works with multiple coronaries - it looks likely that it does, but should make certain.
         //    // Save data necessary to restart the LPN models at coronary surfaces:
         //      // this is done using type-bound public "get" procedures which allow read-only
         //      ! access to the private variables. (we don't actually need these procedures now, as this code is now internal to the multidomain module [it initially was outside, this is copy-paste, so this could be made cleaner!])
         //      open(unit=73,file='coronary_restart.dat',status='old')
         //      read(73,*) !'# ~=Controlled Coronary Model Restart File=~'
         //      do ii = 1,this%numberOfControlledCoronaryModelSurfaces
         //        read(73,*) !'# Coronary outlet index number:', ii
         //        read(73,*) !'# Resistance rp:'
         //        read(73,*) readvalue ! \todo should probably check if the read succeeded so we don't get silent failures here!
         //        call this%setRpByIndex(ii,readvalue)
         //        read(73,*) !'# Resistance rd:'
         //        read(73,*) readvalue
         //        call this%setRdByIndex(ii,readvalue)
         //        read(73,*) !'# Coronary LPN internal pressure P_1:'
         //        read(73,*) readvalue
         //        call this%setP_1ByIndex(ii,readvalue)
         //        read(73,*) !'# Coronary LPN internal pressure P_2:'
         //        read(73,*) readvalue
         //        call this%setP_2ByIndex(ii,readvalue)
         //        read(73,*) !'# currentMyocardialHungerSignal:'
         //        read(73,*) readvalue
         //        call this%setcurrentMeanMyocardialOxygenHungerByIndex(ii,readvalue)
         //        read(73,*) !'# MVO2previousDt:'
         //        read(73,*) readvalue
         //        call this%setMVO2previousByIndex(ii,readvalue)
         //        read(73,*) !'# MVO2'
         //        read(73,*) this%MVO2(ii)
         //        read(73,*) !'# P_a (LPN pressure at outlet itself):'
         //        read(73,*) readvalue
         //        this%P_a(ii) = readvalue
         //        read(73,*) !'# deltaMVO2PerTimestepOverPreviousBeat (used for MVO2 ramping on the present beat)'
         //        read(73,*) readvalue
         //        this%deltaMVO2PerTimestepOverPreviousBeat(ii) = readvalue
         //        read(73,*) !'# Total myocardial oxygen demand during the last complete beat:'
         //        read(73,*) this%MVO2computedAtEndOfPreviousBeat(ii)
         //        read(73,*) !'# Total myocardial oxygen demand during the beat before the last complete beat:'
         //        read(73,*) this%MVO2computedAtEndOfBeatPreviousToPreviousBeat(ii)
         //        read(73,*) !'# Latest myocardial hunger delta ("small h(t)")'
         //        read(73,*) this%hungerDelta(ii)
         //      enddo
         //      read(73,*) !'# Timesteps since the start of the last heart period (for MVO2 calculation)'
         //      read(73,*) hrt%timestepsSinceLastHeartPeriodReset
         //      write(73,*) !'# Are we starting on the first dt of a new heart-beat? 1=yes, 0=no.'
         //      write(73,*) hrt%newBeatJustStarted
         //      close(73)

         //      ! Save O2supplyDemandHistory as a binary file:
         //      ! (easier just to have this done in the multidomain module, and called here)
         //      call this%loadO2supplyDemandHistory()

         //      do ii = 1,this%numberOfControlledCoronaryModelSurfaces
         //        this%oneOverR(ii) = 1d0/(this%rp(ii) + this%rd(ii) + this%ra(ii))
         //      enddo
         // endif
}

void controlledCoronary::computeCapacitorsTopPressures(const double alfi_delt)
{
    // We're now going to solve the 2x2 system m*[P_1;P_2] = rhs.
    // Define the LPN system matrix:
    double m11 = resistanceNearAorta*complianceNearAorta/alfi_delt;
    double m12 = resistanceNearAorta * (1.0/distalResistance + intramyocardialCompliance/alfi_delt);
    double m22 = midResistance * (intramyocardialCompliance/alfi_delt + 1.0/distalResistance) + 1.0;
    //the m21 entry of this matrix is just -1.
    double determinant = m11*m22 + m12;

    double rhs_1 = *flow_n_ptr * resistanceNearAorta + (capacitorNearAortaTopPressure_acceptedAtEndOfLastTimestep*complianceNearAorta + intramyocardialCapacitorTopPressure_acceptedAtEndOfLastTimestep*intramyocardialCompliance + P_IM_mid*intramyocardialCompliance - P_IM_mid_lasttimestep
             *intramyocardialCompliance) * resistanceNearAorta/alfi_delt;
    double rhs_2 = midResistance*intramyocardialCompliance/alfi_delt * (intramyocardialCapacitorTopPressure_acceptedAtEndOfLastTimestep + P_IM_mid - P_IM_mid_lasttimestep);

    // Find the inverse of this matrix
    double inverseM_11 = m22 / determinant;
    double inverseM_12 = -m12 / determinant;
    double inverseM_21 = 1.0 / determinant;
    double inverseM_22 = m11 / determinant;

    // Solve the system for the pressures within the LPN that we need to know:
    capacitorNearAortaTopPressure = inverseM_11 * rhs_1 + inverseM_12 * rhs_2;
    intramyocardialCapacitorTopPressure = inverseM_21 * rhs_1 + inverseM_22 * rhs_2;
}

//    The coronary system-solve in setimplicitcoeff_controlledCoronary relies on knowing the 
//    values of the pressure at various points throughout the LPN model. These need updating
//    each step; this is done by the subroutine updateLPN_coronary.
void controlledCoronary::updateLPN()
{
	double alfi_delt = alfi_local*delt;
	computeCapacitorsTopPressures(alfi_delt);
}

void controlledCoronary::finalizeLPNAtEndOfTimestep()
{
    double alfi_delt = delt;

    computeCapacitorsTopPressures(alfi_delt);

    capacitorNearAortaTopPressure_acceptedAtEndOfLastTimestep = capacitorNearAortaTopPressure;
    intramyocardialCapacitorTopPressure_acceptedAtEndOfLastTimestep = intramyocardialCapacitorTopPressure;
}