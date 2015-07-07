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
	ControlSystemsManager(const double delt)
	: m_delt(delt)
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
	}
	void createParameterController(const parameter_controller_t controllerType, const boost::shared_ptr<NetlistCircuit> boundaryCondition, const int nodeOrComponentIndex);
	void updateAllControlSystems();

	~ControlSystemsManager()
	{
		// Terminate the Python C extensions
		// Py_Finalize();
	}
private:
	std::vector<boost::shared_ptr<AbstractParameterController>> m_controlSystems;
	const double m_delt;
	int m_rank;

};

#endif