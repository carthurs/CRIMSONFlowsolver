#ifndef PARAMETERCONTROLLER_HXX_
#define PARAMETERCONTROLLER_HXX_

#include "CRIMSONPython.hxx"
#include <cmath>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>
#include "timers.hxx"
#include <sstream>
#include <iostream>
#include "datatypesInCpp.hxx"

// This class can be connected to a netlist LPN parameter
// (e.g. resistance for a component, or compliance for a compliance chamber).
// It will adjust that parameter according to necessary control inputs.
class AbstractParameterController
{
public:
	AbstractParameterController(double* const parameterToControl, const int surfaceIndex, const int nodeOrComponentIndex, const circuit_item_t typeOfControlledItem, const int startingTimestepIndex)
	: mp_parameterToControl(parameterToControl), // Set the pointer to the parameter we will control
	m_originalValueOfParameter(*parameterToControl), // Save the original state of the parameter
	m_surfaceIndex(surfaceIndex),
	m_nodeOrComponentIndex(nodeOrComponentIndex),
	m_typeOfControlledItem(typeOfControlledItem),
	m_startingTimestepIndex(startingTimestepIndex)
	{
	}
	virtual ~AbstractParameterController()
	{
	}

	// Function to adjust the controlled parameter as necessary.
	virtual void updateControl() = 0;

	int getIndexOfAssociatedSurface() const;
protected:
	double* const mp_parameterToControl;
	const double m_originalValueOfParameter;
	const int m_surfaceIndex;
	const int m_nodeOrComponentIndex;
	const circuit_item_t m_typeOfControlledItem;
	const int m_startingTimestepIndex;
	virtual void setupControlStateOnRestart() {};
};


class LeftVentricularElastanceController : public AbstractParameterController
{
public:
	LeftVentricularElastanceController(double* const parameterToControl, const int surfaceIndex, const int nodeOrComponentIndex, const circuit_item_t typeOfControlledItem, const double delt, const int startingTimestepIndex)
	: AbstractParameterController(parameterToControl, surfaceIndex, nodeOrComponentIndex, typeOfControlledItem, startingTimestepIndex),
	m_delt(delt)
	{
		m_periodicTime = m_startingTimestepIndex * delt;
		m_timeToMaximumElastance = 0.2782;
		m_timeToRelax = 0.1391;
		m_minimumElastance = 4.10246e-3;
		m_maximumElastance = 3.0827e-1;
		m_heartPeriod = 0.86;
		setupControlStateOnRestart();
	}

	void updateControl() override;
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
	void setupControlStateOnRestart() override;
};

class BleedController : public AbstractParameterController
{
public:
	BleedController(double* const parameterToControl, const int surfaceIndex, const int nodeOrComponentIndex, const circuit_item_t typeOfControlledItem, const int startingTimestepIndex)
	: AbstractParameterController(parameterToControl, surfaceIndex, nodeOrComponentIndex, typeOfControlledItem, startingTimestepIndex)
	{
		int initialTimestep = 0; //\todo sort the restarts for this
		int triggerTimestep = 1600;
		mp_timer = boost::shared_ptr<BasicTimer> (new BasicTimer(initialTimestep, triggerTimestep));
	}
	void updateControl() override;
private:
	boost::shared_ptr<BasicTimer> mp_timer;
	bool m_bleedingOn;
};


class GenericPythonController
{
public:
	GenericPythonController(const double delt, const std::string controllerPythonScriptBaseName, const int startingTimestepIndex)
	: m_delt(PyFloat_FromDouble(delt)),
	m_controllerPythonScriptBaseName(controllerPythonScriptBaseName),
	m_startingTimestepIndex_genericController(startingTimestepIndex)
	{
		// moved to caller because polymorphic function calls dont work in constructors
		// initialise(); 
	}
	virtual ~GenericPythonController()
	{
		safe_Py_DECREF(m_delt);
		safe_Py_DECREF(m_pythonScriptName);
		safe_Py_DECREF(m_pythonControllerClassName);
		safe_Py_DECREF(m_updateControlPyobjectName);
		safe_Py_DECREF(m_customPythonClass);
		safe_Py_DECREF(m_customPythonModule);
		safe_Py_DECREF(m_pythonControllerInstance);
		safe_Py_DECREF(m_pickleMethodName_py);
	}

	void getBroadcastStateData(PyObject*& stateDataBroadcastByThisController);
	void giveStateDataFromOtherPythonControllers(PyObject* allPackagedBroadcastData);
	long getPriority();
	void picklePythonController(const int currentTimestepIndex) const;
	virtual void updateControl();
	void initialise();
protected:

	std::string m_updateControlNameString;
	std::string m_controllerClassName;
	std::string m_passDataToPythonMethodName;
	std::string m_unpickleMethodName;

	PyObject* m_delt;
	PyObject* m_pythonScriptName;
	PyObject* m_pythonControllerClassName;
	PyObject* m_updateControlPyobjectName;
	PyObject* m_customPythonClass;
	PyObject* m_customPythonModule;
	PyObject* m_pythonControllerInstance;
	PyObject* m_pickleMethodName_py;

	const std::string m_controllerPythonScriptBaseName;
	const int m_startingTimestepIndex_genericController;
private:
	virtual std::string getControllerNameQualification();
};


// This class supports user-defined parameter controllers, which the 
// user provides in an external Python script.
class UserDefinedCustomPythonParameterController : public AbstractParameterController, public GenericPythonController
{
public:
	UserDefinedCustomPythonParameterController(double* const parameterToControl, const int surfaceIndex, const int nodeOrComponentIndex, const circuit_item_t typeOfControlledItem, const double delt, const std::string controllerPythonScriptBaseName, const std::vector<std::pair<int,double*>> flowPointerPairs, const std::vector<std::pair<int,double*>> pressurePointerPairs, const std::vector<std::pair<int,double*>> volumePointerPairs, const int startingTimestepIndex)
	: AbstractParameterController(parameterToControl, surfaceIndex, nodeOrComponentIndex, typeOfControlledItem, startingTimestepIndex),
	GenericPythonController(delt, controllerPythonScriptBaseName, startingTimestepIndex),
	m_pressurePointerPairs(pressurePointerPairs),
	m_flowPointerPairs(flowPointerPairs),
	m_volumePointerPairs(volumePointerPairs)
	{
		// moved to caller because polymorphic function calls dont work in constructors
		// initialise();
	}

	void updateControl() override;

	~UserDefinedCustomPythonParameterController() override
	{
		// Py_DECREF(); is like deleting stuff (marking it for deletion by Python garbage
		// collection by decrementing the reference count - deletion happens when
		// the ref count reaches zero, and these calls should make them zero.)
	}

private:
	std::string getControllerNameQualification() override;

	const std::vector<std::pair<int,double*>> m_pressurePointerPairs;
	const std::vector<std::pair<int,double*>> m_flowPointerPairs;
	const std::vector<std::pair<int,double*>> m_volumePointerPairs;

	int errFlag;

};

#endif
