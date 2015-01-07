#ifndef CONTROLLEDCORONARY_HXX_
#define CONTROLLEDCORONARY_HXX_

#include "abstractBoundaryCondition.hxx"
#include "boundaryConditionManager.hxx"
#include "fileReaders.hxx"

class controlledCoronary : public abstractBoundaryCondition
{
public:
	controlledCoronary(int surfaceIndex_in)
	: abstractBoundaryCondition(surfaceIndex_in)
	{
		// Note the index of this RCR (zero-indexed), and count its existance
		// (the order of these two lines is correct!)
		indexOfThisCoronary = numberOfInitialisedCoronaries;
		numberOfInitialisedCoronaries++;
		
		controlledCoronaryReader* controlledCoronaryReader_instance = controlledCoronaryReader::Instance();

		resistanceNearAorta = controlledCoronaryReader_instance->getResistanceNearAorta().at(indexOfThisCoronary);
		complianceNearAorta = controlledCoronaryReader_instance->getComplianceNearAorta().at(indexOfThisCoronary);
		midResistance = controlledCoronaryReader_instance->getMidResistance().at(indexOfThisCoronary);
		intramyocardialCompliance = controlledCoronaryReader_instance->getIntramyocardialCompliance().at(indexOfThisCoronary);
		distalResistance = controlledCoronaryReader_instance->getDistalResistance().at(indexOfThisCoronary);

		minimumAllowedResistance = controlledCoronaryReader_instance->getMinimumAllowedResistance().at(indexOfThisCoronary);
		maximumAllowedResistance = controlledCoronaryReader_instance->getMaximumAllowedResistance().at(indexOfThisCoronary);
		MVO2previousDt = controlledCoronaryReader_instance->getPerfusionBedMVO2_previous().at(indexOfThisCoronary);
		MVO2 = controlledCoronaryReader_instance->getPerfusionBedMVO2_current().at(indexOfThisCoronary);
		proportionOfMyocardiumPerfusedByThisSurface = controlledCoronaryReader_instance->getProportionOfMyocardiumPerfusedByThisSurface().at(indexOfThisCoronary);
		metabolicFeedbackGain = controlledCoronaryReader_instance->getMetabolicFeedbackGain().at(indexOfThisCoronary);
		alphaAdrenergicFeedforwardGain = controlledCoronaryReader_instance->getAlphaAdrenergicFeedforwardGain().at(indexOfThisCoronary);
		betaAdrenergicFeedforwardGain = controlledCoronaryReader_instance->getBetaAdrenergicFeedforwardGain().at(indexOfThisCoronary);
		feedbackDamping = controlledCoronaryReader_instance->getFeedbackDamping().at(indexOfThisCoronary);
		O2supplyDemandHistoryWindowLength_seconds = controlledCoronaryReader_instance->getO2DemandIntegrationWindow().at(indexOfThisCoronary);

		capacitorNearAortaTopPressure_acceptedAtEndOfLastTimestep = controlledCoronaryReader_instance->getCapacitorNearAortaTopPressure().at(indexOfThisCoronary);
		intramyocardialCapacitorTopPressure_acceptedAtEndOfLastTimestep = controlledCoronaryReader_instance->getIntramyocardialCapacitorTopPressure().at(indexOfThisCoronary);

		initialiseModel();

	}

	std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt);

	void updateLPN();

	void finalizeLPNAtEndOfTimestep();

	~controlledCoronary()
	{
		numberOfInitialisedCoronaries--;
	}
private:
	void initialiseModel();
	int indexOfThisCoronary;
	static int numberOfInitialisedCoronaries;

	void computeCapacitorsTopPressures(const double alfi_delt);

	double resistanceNearAorta;
	double complianceNearAorta;
	double midResistance;
	double intramyocardialCompliance;
	double distalResistance;

	double minimumAllowedResistance;
	double maximumAllowedResistance;

	double MVO2previousDt;
	double MVO2;
	double proportionOfMyocardiumPerfusedByThisSurface;
	
	double metabolicFeedbackGain;
	double alphaAdrenergicFeedforwardGain;
	double betaAdrenergicFeedforwardGain;
	double feedbackDamping;
	
	double O2supplyDemandHistoryWindowLength_seconds;
	int O2supplyDemandHistoryWindowLength_timesteps;

	double arterialO2VolumeProportion;
	std::vector<double> O2supplyDemandHistory;
	double currentMyocardialHungerSignal;
	double oneOverTotalResistance;
	double intramyocardialPressureToLVScaling;

	// double inflowPressureInLPN; //\todo why do you store this? Surely you have it in another variable somewhere...
	double capacitorNearAortaTopPressure;
	double intramyocardialCapacitorTopPressure;

	double capacitorNearAortaTopPressure_acceptedAtEndOfLastTimestep;
	double intramyocardialCapacitorTopPressure_acceptedAtEndOfLastTimestep;

	double hungerDelta;

	double P_IM_mid;
	double P_IM_mid_lasttimestep;	

	double deltaMVO2PerTimestepOverPreviousBeat;

	double MVO2computedAtEndOfBeatPreviousToPreviousBeat;
	double MVO2computedAtEndOfPreviousBeat;
};

#endif