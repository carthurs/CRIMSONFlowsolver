#ifndef IMPEDANCEBOUNDARYCONDITION_HXX_
#define IMPEDANCEBOUNDARYCONDITION_HXX_
#include <vector>
#include <deque>
#include <sstream>
#include <string>
#include "AbstractBoundaryCondition.hxx"
#include "datatypesInCpp.hxx"
#include "Interpolators.hxx"

#ifdef BOOST_NO_CXX11_SMART_PTR
#include <boost/shared_ptr.hpp>
#endif

void checkSafetyCounter(const int safetyCounter, const char* currentAction);

class ImpedanceBoundaryCondition : public AbstractBoundaryCondition
{
public:
	ImpedanceBoundaryCondition(const int surfaceIndex_in, const double hstep_in, const double delt_in, const double alfi_in, const double startingTimestepIndex, const int maxsurf, const int nstep)
	: AbstractBoundaryCondition(surfaceIndex_in, hstep_in, delt_in, alfi_in, startingTimestepIndex, maxsurf, nstep),
	m_startingTimestepIndex(startingTimestepIndex)
	{
		{
			std::stringstream filename;
			filename << "Qhistor_impedance_surface_" << surfaceIndex << ".dat";
			m_impedanceFlowHistoryFileName = filename.str();
		}

		{
			std::stringstream filename;
			filename << "time_varying_impedance_surface_" << surfaceIndex << ".dat"; // formerly impt.dat
			m_timeDomainImpedanceFileName = filename.str();
		}

		// if (impfile > 0) // then !for impedance BC #fortran
		// {
		loadInputFiles();
		// }
	}
	void writeImpedanceSpecificFlowHistory() const;
	void resetStateUsingKalmanFilteredEstimate(const double flow, const double pressure, const int timestepNumber) {std::cout << "kalman filter not implemented in ImpedanceBoundaryCondition";};
	void initialiseModel() {}; // not currently used for ImpedanceBoundaryCondition
private:
	const int m_startingTimestepIndex;

	std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt) override;
	void loadInputFiles();
	void readTimeDomainImpedance();
	void readImpedanceFlowHistory();
	void updateStoredFlowHistory();
	std::deque<double> m_qHistImp;
	std::vector<double> m_qHistTry;
	std::vector<double> m_qHistTryF;
	TimeValuePairVector m_timeVaryingImpedance;
	
	#ifdef BOOST_NO_CXX11_SMART_PTR
	boost::shared_ptr<LinearInterpolator> mp_impedanceLinearInterpolator;
	#else
	// prefer unique_ptr if it's available
	std::unique_ptr<LinearInterpolator> mp_impedanceLinearInterpolator;
	#endif

	size_t m_numberOfTimePointsInData;
	std::string m_impedanceFlowHistoryFileName;
	std::string m_timeDomainImpedanceFileName;
};

#endif