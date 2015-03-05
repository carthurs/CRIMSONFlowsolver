#ifndef PUREZERODDRIVER_HXX_
#define PUREZERODDRIVER_HXX_

#include "NetlistBoundaryCondition.hxx"
#include "Netlist3DDomainReplacement.hxx"
#include "boundaryConditionManager.hxx"
#include <boost/shared_ptr.hpp>

class PureZeroDDriver
{
public:
	PureZeroDDriver()
	{
		if (timdat.lstep > 0)
		{
			m_thisIsARestartedSimulation = 1;
	        SimpleFileReader numstartReader("numstart.dat");

	        bool success = false;
	        std::string numstartString = numstartReader.getNextDataSplitBySpacesOrEndOfLine(success);
	        assert(success);

	        m_nextTimestepWrite_zeroDBoundaries_start = boost::lexical_cast<int>(numstartString)+1; // +1 because numstart should contain the step just written before the program last terminated. So we need to start writing on the next (+1 th) time-step.
		}
		else
		{
			m_thisIsARestartedSimulation = 0;
	        m_nextTimestepWrite_zeroDBoundaries_start = 0;
		}
	}
	void init();
	void iter_init();
	void iter_step();
	void iter_finalize();
	void finalize();
private:
	// this is not really a boundary condition here; we just use the machinery of the Netlist to make
	// a replacement for the 3D domain.
	boost::shared_ptr<Netlist3DDomainReplacement> m_zeroDDomainLPN;
	boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
	std::vector<std::pair<boundary_data_t,double>> m_pressuresOrFlowsAtBoundaries;
	double* mp_interfaceFlowsToBeReadByBoundaryConditions;
	double* mp_interfacePressuresToBeReadByBoundaryConditions;

	double* mp_interfaceFlowsToBeReadBy3DDomainReplacement;
	double* mp_interfacePressuresToBeReadBy3DDomainReplacement;

	int m_timestepNumber;
	int m_nextTimestepWrite_zeroDBoundaries_start;
	int m_thisIsARestartedSimulation;

	void placePressuresAndFlowsInStorageArrays_toGiveToBoundaryConditions(std::vector<double> boundaryPressures, std::vector<double> boundaryFlows);
	void placePressuresAndFlowsInStorageArrays_toGiveTo3DDomainReplacement();
};

#endif