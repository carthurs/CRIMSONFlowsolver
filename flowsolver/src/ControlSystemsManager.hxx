#ifndef CONTROLSYSTEMSMANAGER_HXX_
#define CONTROLSYSTEMSMANAGER_HXX_

#include <Python.h>
#include "ParameterController.hxx"
#include "NetlistBoundaryCondition.hxx"
#include "datatypesInCpp.hxx"
#include "CircuitData.hxx"
#include "boost/shared_ptr.hpp"
#include "mpi.h"


// The job of this class is to hold and manage calls to the control systems
// (for example, derived classes of the AbstractParameterController).
// It should be instantiated in boundaryConditionManager, which should
// provide any necessary links to Fortran, and also tell the
// ControlSystemsManager when to do things.
class ControlSystemsManager
{
public:
	ControlSystemsManager(const double delt, const bool masterControlScriptPresent, const int startingTimestepIndex, const int timestepsBetweenRestarts)
	: m_delt(delt),
	m_workingDirectory(boost::filesystem::current_path()),
	m_hasMasterPythonController(masterControlScriptPresent),
	m_startingTimestepIndex(startingTimestepIndex),
	m_timestepsBetweenRestarts(timestepsBetweenRestarts)
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
		if (m_hasMasterPythonController)
		{
			createMasterPythonController();
		}
		m_currentTimestepIndex = m_startingTimestepIndex;
	}
	void createParameterController(const parameter_controller_t controllerType, const boost::shared_ptr<NetlistCircuit> boundaryCondition, const int nodeOrComponentIndex);
	void updateBoundaryConditionControlSystems();
	int getNumberOfControlSystems() const;

	~ControlSystemsManager(){}
private:
	std::vector<boost::shared_ptr<AbstractParameterController>> m_nonPythonControlSystems;
	std::vector<boost::shared_ptr<GenericPythonController>> m_pythonControlSystems; // This vector is for convenience only; it duplicates some elements of m_controlSystems
	std::vector<PyObject*> m_pythonBroadcastDataFromEachController;
	const double m_delt;
	int m_rank;
	void setupPythonBoilerplateScriptPaths();
	void copyFileToWorkingDirectory(const boost::filesystem::path sourcePath, const std::string targetFileName) const;
	boost::filesystem::path m_pathToBoilerplatePythonFlowPrescriberScript;
	boost::filesystem::path m_pathToBoilerplatePythonPressurePrescriberScript;
	const boost::filesystem::path m_workingDirectory;
	void updateAndPassStateInformationBetweenPythonParameterControllers();
	const bool m_hasMasterPythonController;
	boost::shared_ptr<GenericPythonController> mp_masterPythonController;
	void createMasterPythonController();
	void sortPythonControlSystemsByPriority();
	const int m_startingTimestepIndex;
	const int m_timestepsBetweenRestarts;
	int m_currentTimestepIndex;
	void writePythonControlSystemsRestarts();
};

#endif