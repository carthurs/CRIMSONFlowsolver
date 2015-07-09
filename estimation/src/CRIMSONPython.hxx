#ifndef CRIMSONPYTHON_HXX_
#define CRIMSONPYTHON_HXX_

#include "Python.h"

void initialisePython();

void safe_Py_DECREF(PyObject* objectToDecref);

#endif