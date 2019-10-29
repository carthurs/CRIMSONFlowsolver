// Copyright (C) 2008-2010, INRIA
// Author(s): Vivien Mallet, Claire Mouton, Marc Fragu
//
// This file is part of the data assimilation library Verdandi.
//
// Verdandi is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Verdandi. If not, see http://www.gnu.org/licenses/.
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/

%include "verdandi.def"

// SWIG interface to Verdandi.
%module verdandi
%{
#define VERDANDI_DENSE
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#include "VerdandiHeader.hxx"
#include "VerdandiBase.hxx"

#include "model/QuadraticModel.hxx"
#include "model/ClampedBar.hxx"
#include "model/PythonModel.hxx"
#include "observation_manager/GridToNetworkObservationManager.hxx"
#include "observation_manager/LinearObservationManager.hxx"
#include "observation_manager/PythonObservationManager.hxx"
#include "method/OptimalInterpolation.hxx"
#include "method/ForwardDriver.hxx"
#include "method/ReducedOrderExtendedKalmanFilter.hxx"
#include "share/OutputSaver.hxx"
  %}

%include "std_string.i"
%include "std_pair.i"
using namespace std;

%import "seldon/seldon.i"

%rename(SeldonError) Seldon::Error;
%rename(OpsError) Ops::Error;

%include "share/Error.hxx"
%include "share/VerdandiOps.hxx"
%include "seldon/share/Errors.hxx"
%include "ops/Error.hxx"

%exception
{
  try
    {
      $action
	}
  catch(Verdandi::Error& e)
    {
      PyErr_SetString(PyExc_Exception, e.What().c_str());
      return NULL;
    }
  catch(Seldon::Error& e)
    {
      PyErr_SetString(PyExc_Exception, e.What().c_str());
      return NULL;
    }
  catch(Ops::Error& e)
    {
      PyErr_SetString(PyExc_Exception, e.What().c_str());
      return NULL;
    }
  catch(std::exception& e)
    {
      PyErr_SetString(PyExc_Exception, e.what());
      return NULL;
    }
  catch(std::string& s)
    {
      PyErr_SetString(PyExc_Exception, s.c_str());
      return NULL;
    }
  catch(const char* s)
    {
      PyErr_SetString(PyExc_Exception, s);
      return NULL;
    }
  catch(...)
    {
      PyErr_SetString(PyExc_Exception, "Unknown exception...");
      return NULL;
    }
}

%include "VerdandiHeader.hxx"
%include "VerdandiBase.hxx"

%include "model/QuadraticModel.hxx"
%include "model/ClampedBar.hxx"
%include "model/PythonModel.hxx"
%include "observation_manager/GridToNetworkObservationManager.hxx"
%include "observation_manager/LinearObservationManager.hxx"
%include "observation_manager/PythonObservationManager.hxx"
%include "method/OptimalInterpolation.hxx"
%include "method/ForwardDriver.hxx"
%include "method/ReducedOrderExtendedKalmanFilter.hxx"
%include "share/OutputSaver.hxx"

namespace Verdandi
{
  %template(Model) VSWIG_MODEL;
  %template(Model1) VSWIG_MODEL1;
  %template(GTNObservationManager) VSWIG_GRID_TO_NETWORK_OBSERVATION;
  %template(LObservationManager) VSWIG_LINEAR_OBSERVATION;
  %template(Method) VSWIG_METHOD;
  %template(Method1) VSWIG_METHOD1;
  %template(Method2) VSWIG_METHOD2;
  %template(Method3) VSWIG_METHOD3;
  %template(Method4) VSWIG_METHOD4;
  %template(Method5) VSWIG_METHOD5;
}


// For conversions from Seldon to Numpy, and from Numpy to Seldon.
%pythoncode %{
def load_vector(filename, array = False):
    """
    Loads a Seldon vector (in double precision) from a file. If 'array' is set
    to True, the vector is converted to a numpy array.
    """
    import seldon
    if isinstance(filename, str):
        stream = seldon.ifstream(filename)
    else: # assuming 'filename' is already a stream.
        stream = filename

    vector = seldon.VectorDouble()
    vector.Read(stream, True)

    if array:
        import numpy
        vector = numpy.array(vector)

    if isinstance(filename, str):
        stream.close()

    return vector


def load_vector_list(filename, array = False):
    """
    Loads a list of Seldon vectors (in double precision) from a file. If
    'array' is set to True, the vectors are converted to numpy arrays.
    """
    import seldon
    if isinstance(filename, str):
        stream = seldon.ifstream(filename)
    else: # assuming 'filename' is already a stream.
        stream = filename

    vector_list = []
    while stream.peek() != -1:
        vector = seldon.VectorDouble()
        vector.Read(stream, True)
        vector_list.append(vector)

    if array:
        import numpy
        vector_list = [numpy.array(x) for x in vector_list]

    if isinstance(filename, str):
        stream.close()

    return vector_list


def to_vector(v):
    """
    Converts a list or a numpy array to a Seldon vector (in double precision).
    """
    import seldon
    out = seldon.VectorDouble(len(v))
    for i in range(len(v)):
        out[i] = v[i]
    return out


def load_matrix(filename, array = False):
    """
    Loads a Seldon matrix (in double precision) from a file. If 'array' is set
    to True, the matrix is converted to a numpy array.
    """
    import seldon
    if isinstance(filename, str):
        stream = seldon.ifstream(filename)
    else: # assuming 'filename' is already a stream.
        stream = filename

    matrix = seldon.MatrixDouble()
    matrix.Read(stream, True)

    if array:
        import numpy
        matrix = numpy.array(matrix)

    if isinstance(filename, str):
        stream.close()

    return matrix


def load_matrix_list(filename, array = False):
    """
    Loads a list of Seldon matrices (in double precision) from a file. If
    'array' is set to True, the matrices are converted to numpy arrays.
    """
    import seldon
    if isinstance(filename, str):
        stream = seldon.ifstream(filename)
    else: # assuming 'filename' is already a stream.
        stream = filename

    matrix_list = []
    while stream.peek() != -1:
        matrix = seldon.MatrixDouble()
        matrix.Read(stream, True)
        matrix_list.append(matrix)

    if array:
        import numpy
        matrix_list = [numpy.array(x) for x in matrix_list]

    if isinstance(filename, str):
        stream.close()

    return matrix_list


def to_matrix(m):
    """
    Converts a numpy array to a Seldon matrix (in double precision).
    """
    import seldon
    out = seldon.MatrixDouble(m.shape[0], m.shape[1])
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            out[i, j] = m[i, j]
    return out
%}
