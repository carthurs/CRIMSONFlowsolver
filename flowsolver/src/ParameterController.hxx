#ifndef PARAMETERCONTROLLER_HXX_
#define PARAMETERCONTROLLER_HXX_

#include <Python.h>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>
#include "timers.hxx"
#include <sstream>
#include <iostream>

// This class can be connected to a netlist LPN parameter
// (e.g. resistance for a component, or compliance for a compliance chamber).
// It will adjust that parameter according to necessary control inputs.
class AbstractParameterController
{
public:
	AbstractParameterController(double* const parameterToControl)
	: mp_parameterToControl(parameterToControl), // Set the pointer to the parameter we will control
	m_originalValueOfParameter(*parameterToControl) // Save the original state of the parameter
	{
	}
	virtual ~AbstractParameterController()
	{
	}

	// Function to adjust the controlled parameter as necessary.
	virtual void updateControl() = 0;
protected:
	double* const mp_parameterToControl;
	const double m_originalValueOfParameter;
};


class LeftVentricularElastanceController : public AbstractParameterController
{
public:
	LeftVentricularElastanceController(double* const parameterToControl, const double delt)
	: AbstractParameterController(parameterToControl),
	m_delt(delt)
	{
		m_periodicTime = 0.0; //\todo think about this for restarts!
		m_timeToMaximumElastance = 0.2782;
		m_timeToRelax = 0.1391;
		m_minimumElastance = 4.10246e-3;
		m_maximumElastance = 3.0827e-1;
		m_heartPeriod = 0.86;
	}

	void updateControl();
private:
	const double m_delt;
	
	double m_periodicTime;
	double m_heartPeriod;
	
	double m_timeToMaximumElastance;
	double m_timeToRelax;

	double m_minimumElastance;
	double m_maximumElastance;

	double getElastance();

	void updatePeriodicTime();
};

class BleedController : public AbstractParameterController
{
public:
	BleedController(double* const parameterToControl)
	: AbstractParameterController(parameterToControl)
	{
		int initialTimestep = 0; //\todo sort the restarts for this
		int triggerTimestep = 1600;
		mp_timer = boost::shared_ptr<BasicTimer> (new BasicTimer(initialTimestep, triggerTimestep));
	}
	void updateControl();
private:
	boost::shared_ptr<BasicTimer> mp_timer;
	bool m_bleedingOn;
};

// This class supports user-defined parameter controllers, which the 
// user provides in an external Python script.
//
// Currently, it only supports
class UserDefinedCustomPythonParameterController : public AbstractParameterController
{
public:
	UserDefinedCustomPythonParameterController(double* const parameterToControl, const double delt, const std::string controllerPythonScriptBaseName, const std::vector<std::pair<int,double*>> flowPointerPairs, const std::vector<std::pair<int,double*>> pressurePointerPairs, const std::vector<std::pair<int,double*>> volumePointerPairs)
	: AbstractParameterController(parameterToControl),
	m_delt(PyFloat_FromDouble(delt)),
	m_controllerPythonScriptBaseName(controllerPythonScriptBaseName),
	m_pressurePointerPairs(pressurePointerPairs),
	m_flowPointerPairs(flowPointerPairs),
	m_volumePointerPairs(volumePointerPairs)
	{
		initialise();
	}


	~UserDefinedCustomPythonParameterController()
	{
		// Py_DECREF(); is like deleting stuff (marking it for deletion by Python garbage
		// collection by decrementing the reference count - deletion happens when
		// the ref count reaches zero, and these calls should make them zero.)
		safe_Py_DECREF(m_delt);
		safe_Py_DECREF(m_pythonScriptName);
		safe_Py_DECREF(m_pythonControllerClassName);
		safe_Py_DECREF(m_updateControlPyobjectName);
		safe_Py_DECREF(m_customPythonClass);
		safe_Py_DECREF(m_customPythonModule);
		safe_Py_DECREF(m_pythonParameterControllerInstance);
	}

	void updateControl();
private:
	PyObject* m_delt;
	PyObject* m_pythonScriptName;
	PyObject* m_pythonControllerClassName;
	PyObject* m_updateControlPyobjectName;
	PyObject* m_customPythonClass;
	PyObject* m_customPythonModule;
	PyObject* m_pythonParameterControllerInstance;

	std::string m_updateControlNameString;
	std::string m_controllerPythonScriptBaseName;

	const std::vector<std::pair<int,double*>> m_pressurePointerPairs;
	const std::vector<std::pair<int,double*>> m_flowPointerPairs;
	const std::vector<std::pair<int,double*>> m_volumePointerPairs;

	int errFlag;

	void initialise();
	void safe_Py_DECREF(PyObject* toBeDeleted);

};

#endif