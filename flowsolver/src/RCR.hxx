#ifndef RCR_HXX_
#define RCR_HXX_

#include "AbstractBoundaryCondition.hxx"
#include "fileReaders.hxx"
#include "../../estimation/src/SimvascularGlobalArrayTransfer.h"

class RCR : public AbstractBoundaryCondition
{
public:
	RCR(const int surfaceIndex_in, const double hstep_in, const double delt_in, const double alfi_in, const double currentTimestepIndex, const int maxsurf, const int nstep)
	: AbstractBoundaryCondition(surfaceIndex_in, hstep_in, delt_in, alfi_in, currentTimestepIndex, maxsurf, nstep)
	{
		// Note the index of this RCR (zero-indexed), and count its existance
		// (the order of these two lines is correct!)
		indexOfThisRCR = s_numberOfInitialisedRCRs;
		s_numberOfInitialisedRCRs++;
		m_needsPressureToBeInitialisedFromFortran = false;
		
		initialiseModel();
		RcrtReader* rcrtReader_instance = RcrtReader::Instance();
		proximalResistance = rcrtReader_instance->getR1()[indexOfThisRCR];
		capacitance = rcrtReader_instance->getC()[indexOfThisRCR];
		distalResistance = rcrtReader_instance->getR2()[indexOfThisRCR];
		timeDataPdist = rcrtReader_instance->getTimeDataPdist()[indexOfThisRCR];
		lengthOftimeDataPdist = rcrtReader_instance->getNumDataRCR()[indexOfThisRCR];

		// Set up for Kalman filtering:
		SimvascularGlobalArrayTransfer::Get()->setPointerToWindkesselProximalResistance(&proximalResistance, indexOfThisRCR);
		SimvascularGlobalArrayTransfer::Get()->setPointerToWindkesselDistalResistance(&distalResistance, indexOfThisRCR);
		SimvascularGlobalArrayTransfer::Get()->setPointerToWindkesselCompilance(&capacitance, indexOfThisRCR);
		SimvascularGlobalArrayTransfer::Get()->setPointerToRCRSurfacePressure(&pressure_n, indexOfThisRCR);
	}
	
 	std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timeAtStepNplus1, const double alfi_delt) override;

 	void setPressureFromFortran();
 	void resetStateUsingKalmanFilteredEstimate(const double flow, const double pressure, const int timestepNumber) override;
 	void getPressureAndFlowPointersFromFortran() override;


	
//  	procedure :: setimplicitcoeff_rcr => setimplicitcoeff_rcr
        // void setImplicitCoeff();
//     procedure :: updxvars_rcr => updxvars_rcr         
//     procedure :: writexvars_rcr => writexvars_rcr
//     procedure :: assign_ptrs_ext_rcr => assign_ptrs_ext_rcr

	~RCR() override
	{
		s_numberOfInitialisedRCRs--;
	}
protected:
	double linInterpolateTimeData(const double &currentTime, const int timeDataLength);

private:
	void initialiseModel() override;
	static int s_numberOfInitialisedRCRs;
	int indexOfThisRCR;
	int pdmax;
	double proximalResistance; // Proximal resistance
	double capacitance; // Capacitance (compliance)
	double distalResistance; // Distal Resistance
	std::vector<std::pair<double,double>> timeDataPdist; // Time-varying disal pressure data
	int lengthOftimeDataPdist;
	bool m_needsPressureToBeInitialisedFromFortran;
};

#endif