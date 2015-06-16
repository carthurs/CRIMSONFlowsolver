#include "ParameterController.hxx"
#include <iostream>

void LeftVentricularElastanceController::updateControl()
{
	updatePeriodicTime();
	// adjust the controlled elastance:
	*mp_parameterToControl = getElastance();
}

double LeftVentricularElastanceController::getElastance()
{
	// *** analytical elastance function from:
	//     pope, s. r.; ellwein, l. m.; zapata, c. l.; novak, v.; kelley, c. t. & olufsen, m. s.  
	//     estimation and identification of parameters in a lumped cerebrovascular model.
	//     math biosci eng, 2009, 6, 93-115
	
	double elastance;

	// This is the elastance function. It's defined piecewise:
	if ( m_periodicTime <= m_timeToMaximumElastance )
	{
 		elastance = m_minimumElastance
            + 0.5*(m_maximumElastance - m_minimumElastance)
            * (1.0 - cos((m_periodicTime*M_PI)/m_timeToMaximumElastance));
    }
  	else if ( m_periodicTime <= (m_timeToMaximumElastance + m_timeToRelax) )
  	{
     	elastance = m_minimumElastance
            + 0.5*(m_maximumElastance-m_minimumElastance)
            * (1.0 + cos((m_periodicTime-m_timeToMaximumElastance)*(M_PI/m_timeToRelax)));
    }
	else if ( m_periodicTime > (m_timeToMaximumElastance + m_timeToRelax) )
	{
		elastance = m_minimumElastance;
	}

	// std::cout << "elastance was: "<< elastance << std::endl;

  	return elastance;
}

void LeftVentricularElastanceController::updatePeriodicTime()
{
	m_periodicTime = m_periodicTime + m_delt;
	// Keep m_periodicTime in the range [0,m_heartPeriod) :
	if (m_periodicTime >= m_heartPeriod)
	{
		m_periodicTime = m_periodicTime - m_heartPeriod;
	}
	// std::cout << "m_periodicTime was: "<< m_periodicTime << std::endl;
}

void BleedController::updateControl()
{
	bool m_bleedingOn = mp_timer->hasTheTimeCome();
	if (m_bleedingOn)
	{
		*mp_parameterToControl = 0.001; // set the resistance / compliance to be tiny (depending on the type of component we're controlling here...)
	}
	mp_timer->incrementTimer();
}

void UserDefinedCustomPythonParameterController::initialise()
{
	// Catch and report any Python errors:
	try
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
		PyObject* sysPath = PySys_GetObject((char*)"path"); // this is a BORROWED reference, so do not DECREF it!
                PyList_Append(sysPath, PyString_FromString(current_path));
		int success = PySys_SetObject("path",  sysPath);
		if (success != 0)
 		{
			throw std::runtime_error("EE: Failed to append the working directory to the Python path. Contact the developers.\n");
		}

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
	catch(...) // Catch any exception
	{
		if (PyErr_Occurred())
		{
			std::cout << std::endl << "EE: An error occurred when constructing the Python parameter controller "  << std::endl
					               << m_controllerPythonScriptBaseName << ".py. Details are below:, but if this"  << std::endl
					               << " is an error with a missing module, you may need to set paths in your "  << std::endl
					               << "Python script with sys.path.append(\"/path/to/missing/module/\")." << std::endl << std::endl;
			PyErr_Print();
		}
		// Rethrow the original exception (whether or not it was Python's).
		throw;
	}
}

void UserDefinedCustomPythonParameterController::updateControl()
{
	// Catch and report any Python errors:
	try
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
	catch(...) // Catch any exception
	{
		if (PyErr_Occurred())
		{
			std::cout << std::endl << "EE: An error occurred when constructing the Python parameter controller "
					  << m_controllerPythonScriptBaseName << ".py. Details below:" << std::endl;
			PyErr_Print();
		}
		// Rethrow the original exception (whether or not it was Python's).
		throw;
	}
}

void UserDefinedCustomPythonParameterController::safe_Py_DECREF(PyObject* toBeDeleted)
{
	// Avoid trying to delete null pointers:
	if (toBeDeleted)
	{
		Py_DECREF(toBeDeleted);
		toBeDeleted = NULL;
	}
}
