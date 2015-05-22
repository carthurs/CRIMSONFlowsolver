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
		m_updateControlNameString = "updateControl";

		std::stringstream fullFileName;
		fullFileName << m_controllerPythonScriptBaseName << ".py";
		boost::filesystem::path fullFileName_path(fullFileName.str());
		if (!boost::filesystem::exists(fullFileName_path))
		{
			std::stringstream errorMessage;
			errorMessage << "EE: Could not find custom parameter control script " << m_controllerPythonScriptBaseName.c_str() << ".py" << std::endl;
			throw std::runtime_error(errorMessage.str());
		}

		// Change Python's current path to be the same as that which C++ is
		// currently using:
		boost::filesystem::path currentDirectory( boost::filesystem::current_path() );
		char* current_path = (char*) currentDirectory.string().c_str();
		PySys_SetPath(current_path);

		// This is the name of the method that gets called on the class to update the control
		// on each time-step
		m_updateControlPyobjectName = PyString_FromString(m_updateControlNameString.c_str());
		m_pythonControllerClassName = PyString_FromString(m_controllerPythonScriptBaseName.c_str());
		m_pythonScriptName = PyString_FromString(m_controllerPythonScriptBaseName.c_str());
		m_customPythonModule = PyImport_Import(m_pythonScriptName);

		if (m_customPythonModule == NULL)
		{
			std::stringstream errorMessage;
			errorMessage << "EE: Error while parsing file ";
			errorMessage << m_controllerPythonScriptBaseName.c_str() << ".py" << std::endl;
			throw std::runtime_error(errorMessage.str());
		}

		// Get a reference to the custom controller class from within the user-provide Python script
		m_customPythonClass = PyObject_GetAttr(m_customPythonModule, m_pythonControllerClassName);
		if (m_customPythonClass == NULL)
		{
			std::stringstream errorMessage;
			errorMessage << "EE: Could not find a class named " << m_controllerPythonScriptBaseName.c_str();
			errorMessage << " in file " << m_controllerPythonScriptBaseName.c_str() << ".py" << std::endl;
			throw std::runtime_error(errorMessage.str());
		}


		// Instantiate the Python controller class:
		if (PyCallable_Check(m_customPythonClass) == 1)
		{
			// Instantiate the controller class
			m_pythonParameterControllerInstance = PyObject_CallObject(m_customPythonClass, NULL);
		}
		else
		{
			std::stringstream errorMessage;
			errorMessage << "EE: Failed to call a class named " << m_controllerPythonScriptBaseName.c_str();
			errorMessage << " in file " << m_controllerPythonScriptBaseName.c_str() << ".py" << std::endl;
			throw std::runtime_error(errorMessage.str());
		}

	}

	void safe_Py_DECREF(PyObject* toBeDeleted)
	{
		// Avoid trying to delete null pointers:
		if (toBeDeleted)
		{
			Py_DECREF(toBeDeleted);
			toBeDeleted = NULL;
		}
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

	void updateControl()
	{
		// Pack up the pressures and flows in Python dictionaries for this Netlist,
		// indexed by the input data indices for the nodes / componnents:
		PyObject* pressuresInThisNetlist = PyDict_New();
		for (auto pressurePair = m_pressurePointerPairs.begin(); pressurePair != m_pressurePointerPairs.end(); pressurePair++)
		{
			PyObject* nodeIndexInInputData = PyInt_FromLong((long) pressurePair->first);
			PyObject* pressurePointer = PyFloat_FromDouble(*(pressurePair->second));
			errFlag = PyDict_SetItem(pressuresInThisNetlist, nodeIndexInInputData, pressurePointer);
			assert(errFlag == 0);

			safe_Py_DECREF(nodeIndexInInputData);
			safe_Py_DECREF(pressurePointer);
		}

		PyObject* flowsInThisNetlist = PyDict_New();
		for (auto flowPair = m_flowPointerPairs.begin(); flowPair != m_flowPointerPairs.end(); flowPair++)
		{
			PyObject* componentIndexInInputData = PyInt_FromLong((long) flowPair->first);
			PyObject* flowPointer = PyFloat_FromDouble(*(flowPair->second));
			errFlag = PyDict_SetItem(flowsInThisNetlist, componentIndexInInputData, flowPointer);
			assert(errFlag == 0);

			safe_Py_DECREF(componentIndexInInputData);
			safe_Py_DECREF(flowPointer);
		}

		PyObject* volumesInThisNetlist = PyDict_New();
		for (auto volumePair = m_volumePointerPairs.begin(); volumePair != m_volumePointerPairs.end(); volumePair++)
		{
			PyObject* componentIndexInInputData = PyInt_FromLong((long) volumePair->first);
			PyObject* volumePointer = PyFloat_FromDouble(*(volumePair->second));
			errFlag = PyDict_SetItem(volumesInThisNetlist, componentIndexInInputData, volumePointer);
			assert(errFlag == 0);

			safe_Py_DECREF(componentIndexInInputData);
			safe_Py_DECREF(volumePointer);
		}


		// Convert the parameter value to Python format, for passing to Python:
		PyObject* parameterValue = PyFloat_FromDouble(*mp_parameterToControl);
		// Call the updateControl method in the Python script:
		PyObject* newParameterValue = PyObject_CallMethodObjArgs(m_pythonParameterControllerInstance, m_updateControlPyobjectName, parameterValue, m_delt, pressuresInThisNetlist, flowsInThisNetlist, volumesInThisNetlist, NULL);
		if (newParameterValue == NULL)
		{
			std::stringstream errorMessage;
			errorMessage << "EE: Failed to call a method named " << m_updateControlNameString;
			errorMessage << " of a class named " << m_controllerPythonScriptBaseName.c_str();
			errorMessage << " in file " << m_controllerPythonScriptBaseName.c_str() << ".py" << std::endl;
			throw std::runtime_error(errorMessage.str());
		}
		// Place the newly-computed parameter value object being controlled:
		*mp_parameterToControl = PyFloat_AsDouble(newParameterValue);

		safe_Py_DECREF(parameterValue);
		safe_Py_DECREF(newParameterValue);
		safe_Py_DECREF(pressuresInThisNetlist);
		safe_Py_DECREF(flowsInThisNetlist);
		safe_Py_DECREF(volumesInThisNetlist);
	}
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
};

#endif