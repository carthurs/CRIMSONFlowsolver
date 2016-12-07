#ifndef RCR_HXX_
#define RCR_HXX_

#include "AbstractBoundaryCondition.hxx"
#include "fileReaders.hxx"
#include "../../estimation/src/SimvascularGlobalArrayTransfer.h"
#include "datatypesInCpp.hxx"
#include "Interpolators.hxx"

class RCR : public AbstractBoundaryCondition
{
	FRIEND_TEST(testMultidom, checkRCRLinearInterpolators);
public:
	RCR(const int surfaceIndex_in, const double hstep_in, const double delt_in, const double alfi_in, const double currentTimestepIndex, const int maxsurf, const int nstep)
	: AbstractBoundaryCondition(surfaceIndex_in, hstep_in, delt_in, alfi_in, currentTimestepIndex, maxsurf, nstep),
	m_indexOfThisRCR(s_numberOfInitialisedRCRs),
	m_timeDataPdist(RcrtReader::Instance()->getTimeDataPdist()[m_indexOfThisRCR]),
	m_pDistLinearInterpolator(m_timeDataPdist)
	{
		// Note the index of this RCR (zero-indexed), and count its existance
		// (the order of these two lines is correct!)
		s_numberOfInitialisedRCRs++;
		m_needsPressureToBeInitialisedFromFortran = false;
		
		initialiseModel();
		RcrtReader* rcrtReader_instance = RcrtReader::Instance();
		proximalResistance = rcrtReader_instance->getR1()[m_indexOfThisRCR];
		capacitance = rcrtReader_instance->getC()[m_indexOfThisRCR];
		distalResistance = rcrtReader_instance->getR2()[m_indexOfThisRCR];
		lengthOftimeDataPdist = rcrtReader_instance->getNumDataRCR()[m_indexOfThisRCR];

		// Set up for Kalman filtering:
		SimvascularGlobalArrayTransfer::Get()->setPointerToWindkesselProximalResistance(&proximalResistance, m_indexOfThisRCR);
		SimvascularGlobalArrayTransfer::Get()->setPointerToWindkesselDistalResistance(&distalResistance, m_indexOfThisRCR);
		SimvascularGlobalArrayTransfer::Get()->setPointerToWindkesselCompilance(&capacitance, m_indexOfThisRCR);
		SimvascularGlobalArrayTransfer::Get()->setPointerToRCRSurfacePressure(&pressure_n, m_indexOfThisRCR);
	}
	
 	void setPressureFromFortran();
 	void resetStateUsingKalmanFilteredEstimate(const double flow, const double pressure, const int timestepNumber) override;
 	void getPressureAndFlowPointersFromFortran() override;


	void finaliseAtEndOfTimestep() override { throw std::logic_error("Call to unused method RCR::finaliseAtEndOfTimestep(). Contact the developers.");}
//  	procedure :: setimplicitcoeff_rcr => setimplicitcoeff_rcr
        // void setImplicitCoeff();
//     procedure :: updxvars_rcr => updxvars_rcr         
//     procedure :: writexvars_rcr => writexvars_rcr
//     procedure :: assign_ptrs_ext_rcr => assign_ptrs_ext_rcr

	~RCR() override
	{
		s_numberOfInitialisedRCRs--;
	}
private:
	void initialiseModel() override;
	static int s_numberOfInitialisedRCRs;
	
	const int m_indexOfThisRCR;
	const TimeValuePairVector m_timeDataPdist; // Time-varying disal pressure data
	const LinearInterpolator m_pDistLinearInterpolator;

	int pdmax;
	double proximalResistance; // Proximal resistance
	double capacitance; // Capacitance (compliance)
	double distalResistance; // Distal Resistance
	int lengthOftimeDataPdist;
	bool m_needsPressureToBeInitialisedFromFortran;
	std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timeAtStepNplus1, const double alfi_delt) override;
};

#endif