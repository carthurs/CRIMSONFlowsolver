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
	ControlSystemsManager(const double delt, const bool masterControlScriptPresent)
	: m_delt(delt),
	m_workingDirectory(boost::filesystem::current_path()),
	m_hasMasterPythonController(masterControlScriptPresent)
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
		if (m_hasMasterPythonController)
		{
			createMasterPythonController();
		}
	}
	void createParameterController(const parameter_controller_t controllerType, const boost::shared_ptr<NetlistCircuit> boundaryCondition, const int nodeOrComponentIndex);
	void updateBoundaryConditionControlSystems();

	~ControlSystemsManager()
	{
		// Terminate the Python C extensions
		// Py_Finalize();
	}
private:
	std::vector<boost::shared_ptr<AbstractParameterController>> m_controlSystems;
	std::vector<boost::shared_ptr<UserDefinedCustomPythonParameterController>> m_pythonControlSystems; // This vector is for convenience only; it duplicates some elements of m_controlSystems
	std::vector<PyObject*> m_pythonBroadcastDataFromEachController;
	const double m_delt;
	int m_rank;
	void setupPythonBoilerplateScriptPaths();
	void copyFileToWorkingDirectory(const boost::filesystem::path sourcePath, const std::string targetFileName) const;
	boost::filesystem::path m_pathToBoilerplatePythonFlowPrescriberScript;
	boost::filesystem::path m_pathToBoilerplatePythonPressurePrescriberScript;
	const boost::filesystem::path m_workingDirectory;
	void updateAndPassStateInformationBetweenPythonParameterControllers();
	bool m_hasMasterPythonController;
	boost::shared_ptr<GenericPythonController> mp_masterPythonController;
	void createMasterPythonController();

};

#endif