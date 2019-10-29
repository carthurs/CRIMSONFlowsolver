// Copyright (C) 2011-2012 INRIA
// Author(s): Kévin Charpentier, Vivien Mallet
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


#ifndef VERDANDI_FILE_OBSERVATION_MANAGER_PYTHONOBSERVATIONMANAGER_CXX


#include "PythonObservationManager.hxx"


namespace Verdandi
{


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////


    //! Constructor.
    PythonObservationManager::PythonObservationManager()
    {
        Py_Initialize();
        import_array();
        is_module_initialized_ = false;
    }


    //! Destructor.
    PythonObservationManager::~PythonObservationManager()
    {
        Py_Finalize();
    }


    ////////////////////
    // INITIALIZATION //
    ////////////////////


    //! Initializes the observation manager.
    /*!
      \param[in] model model.
      \param[in] configuration_file configuration file.
      \tparam Model the model type; e.g. ShallowWater<double>
    */
    template <class Model>
    void PythonObservationManager::Initialize(const Model& model,
                                              string configuration_file)
    {
        VerdandiOps configuration(configuration_file);

        configuration.Set("python_observation_manager.module", module_);
        configuration.Set("python_observation_manager.directory", directory_);
        configuration.Set("python_observation_manager.class_name",
                          class_name_);

        pyObservationManagerFile_
            = PyString_FromString(module_.c_str());

        if (directory_ != "")
        {
            string PythonPath = "sys.path.append(\"" + directory_ + "\")";
            PyRun_SimpleString("import sys");
            PyRun_SimpleString(PythonPath.c_str());
        }

        pyObservationManagerModule_
            = PyImport_Import(pyObservationManagerFile_);

        pyObservationManagerDict_
            = PyModule_GetDict(pyObservationManagerModule_);

        pyObservationManagerClass_
            = PyDict_GetItemString(pyObservationManagerDict_,
                                   class_name_.c_str());

        PyObject *file_path
            = PyString_FromString(configuration.GetFilePath().c_str());
        PyObject *args = PyTuple_New(2);
        PyTuple_SetItem(args, 0, file_path);
        PyTuple_SetItem(args, 1, PyInt_FromLong(model.GetNstate()));

        if (!PyCallable_Check(pyObservationManagerClass_))
            throw ErrorConfiguration("PythonObservationManager"
                                     "::PythonObservationManager",
                                     "The class \"" + class_name_
                                     + "\" is not defined in \"" + module_
                                     + "\".");

        pyObservationManagerInstance_
            = PyObject_CallObject(pyObservationManagerClass_, args);

        is_module_initialized_ = true;

        Py_DECREF(pyObservationManagerFile_);
        Py_DECREF(pyObservationManagerModule_);
    }


    //! Activates or deactivates the option 'discard_observation'.
    /*!
      \param[in] discard_observation if set to true, each observation will be
      used at most one time.
    */
    void PythonObservationManager
    ::DiscardObservation(bool discard_observation)
    {
        char function_name[] = "DiscardObservation";
        char format_unit[] = "i";

        if (PyObject_CallMethod(pyObservationManagerInstance_,
                                function_name, format_unit,
                                discard_observation) != Py_None)
            throw ErrorPythonUndefined("PythonObservationManager"
                                       "::DiscardObservation",
                                       string(function_name), "(self, time)",
                                       module_);
    }


    //! Sets the time of observations to be loaded.
    /*!
      \param[in] model the model.
      \param[in] time a given time.
    */
    template <class Model>
    void PythonObservationManager
    ::SetTime(const Model& model, double time)
    {
        char function_name[] = "SetTime";
        char format_unit[] = "d";
        if (PyObject_CallMethod(pyObservationManagerInstance_, function_name,
                                format_unit, time) != Py_None)
            throw ErrorPythonUndefined("PythonObservationManager::SetTime",
                                       string(function_name), "(self, time)",
                                       module_);
    }


    /////////////////
    // OBSERVATION //
    /////////////////


    //! Returns the observations.
    /*! This method is called after 'SetTime' set the time at which the
      observations are requested.
      \param[out] observation observation vector.
    */
    void PythonObservationManager
    ::GetObservation(PythonObservationManager::observation& observation)
    {
        char function_name[] = "GetObservation";

        PyObject *pyObs = PyObject_CallMethod(pyObservationManagerInstance_,
                                              function_name, NULL);

        if (pyObs == NULL)
            throw ErrorPythonUndefined("PythonModel::GetObservation",
                                       string(function_name), "(self)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyObs))
            throw ErrorProcessing("PythonModel::GetObservation",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *obsArray = reinterpret_cast<PyArrayObject*>(pyObs);

        int size = obsArray->dimensions[0];
        observation.Reallocate(size);
        memcpy(observation.GetDataVoid(), obsArray->data,
               size * sizeof(double));
    }


    ////////////////
    // INNOVATION //
    ////////////////


    //! Returns an innovation.
    /*! This method is called after 'SetTime' set the time at which the
      innovation is requested.
      \param[in] state state vector.
      \param[out] innovation innovation vector.
    */
    template <class state>
    void PythonObservationManager
    ::GetInnovation(const state& x,
                    PythonObservationManager::observation& innovation)
    {
        char function_name[] = "GetInnovation";
        char format_unit[] = "O";
        npy_intp dim[1];
        dim[0] = x.GetLength();

        PyObject *pyState
            = PyArray_SimpleNewFromData(1, dim, NPY_DOUBLE,
                                        x.GetDataVoid());

        PyObject *pyInnovation
            = PyObject_CallMethod(pyObservationManagerInstance_,
                                  function_name, format_unit, pyState);

        if (pyInnovation == NULL)
            throw ErrorPythonUndefined("PythonModel::GetInnovation",
                                       string(function_name), "(self, state)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyInnovation))
            throw ErrorProcessing("PythonModel::GetObservation",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *innovationArray =
            reinterpret_cast<PyArrayObject*>(pyInnovation);

        int size = innovationArray->dimensions[0];
        innovation.Reallocate(size);
        memcpy(innovation.GetDataVoid(), innovationArray->data,
               size * sizeof(double));
    }


    ////////////
    // ACCESS //
    ////////////


    //! Indicates if some observations are available at a given time.
    /*!
      \param[in] time a given time.
    */
    bool PythonObservationManager::HasObservation(double time)
    {
        char function_name[] = "HasObservation";
        char format_unit[] = "d";
        PyObject *has_observation =
            PyObject_CallMethod(pyObservationManagerInstance_,
                                function_name, format_unit, time);

        if (has_observation == NULL)
            throw ErrorPythonUndefined("PythonModel::HasObservation",
                                       string(function_name), "(self, time)",
                                       module_);

        return PyInt_AsLong(has_observation) != 0;
    }


    //! Indicates if some observations are available at current time.
    bool PythonObservationManager::HasObservation() const
    {
        char function_name[] = "HasObservation";
        PyObject *has_observation =
            PyObject_CallMethod(pyObservationManagerInstance_,
                                function_name, NULL);

        if (has_observation == NULL)
            throw ErrorPythonUndefined("PythonModel::HasObservation",
                                       string(function_name), "(self)",
                                       module_);

        return PyInt_AsLong(has_observation) != 0;
    }


    //! Returns the number of available observations.
    /*!
      \return The total number of observation at current time.
    */
    int PythonObservationManager::GetNobservation() const
    {
        char function_name[] = "GetNobservation";
        PyObject *Nobservation =
            PyObject_CallMethod(pyObservationManagerInstance_,
                                function_name, NULL);

        if (Nobservation == NULL)
            throw ErrorPythonUndefined("PythonModel::GetNobservation",
                                       string(function_name), "(self)",
                                       module_);

        return PyInt_AsLong(Nobservation);
    }


    ///////////////
    // OPERATORS //
    ///////////////


    //! Applies the observation operator to a given vector.
    /*! This method is called after 'SetTime' set the time at which the
      operator is defined.
      \param[in] x a vector.
      \param[out] y the value of the operator applied to \a x. It is resized
      if needed.
    */
    template <class state>
    void PythonObservationManager
    ::ApplyOperator(const state& x, observation& y) const
    {
        char function_name[] = "ApplyOperator";
        char format_unit[] = "O";
        npy_intp dim[1];
        dim[0] = x.GetLength();

        PyObject *pyState =  PyArray_SimpleNewFromData(1, dim, NPY_DOUBLE,
                                                       x.GetDataVoid());

        PyObject *pyY = PyObject_CallMethod(pyObservationManagerInstance_,
                                            function_name,
                                            format_unit, pyState);

        if (pyY == NULL)
            throw ErrorPythonUndefined("PythonObservationManager"
                                       "::ApplyOperator",
                                       string(function_name),
                                       "(self, state, observation)", module_);

        if (!PyArray_ISCONTIGUOUS(pyY))
            throw ErrorProcessing("PythonObservationManager::ApplyOperator",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *yArray = reinterpret_cast<PyArrayObject*>(pyY);

        int size = yArray->dimensions[0];
        y.Reallocate(size);
        memcpy(y.GetDataVoid(), yArray->data, size * sizeof(double));
    }


    //! Applies the tangent linear operator to a given vector.
    /*! This method is called after 'SetTime' set the time at which the
      operator is defined.
      \param[in] x a vector.
      \param[out] y the value of the tangent linear operator applied to \a
      x. It is resized if needed.
    */
    template <class state>
    void PythonObservationManager
    ::ApplyTangentLinearOperator(const state& x, observation& y) const
    {
        char function_name[] = "ApplyTangentLinearOperator";
        char format_unit[] = "O";
        npy_intp dim[1];
        dim[0] = x.GetLength();

        PyObject *pyState =  PyArray_SimpleNewFromData(1, dim, NPY_DOUBLE,
                                                       x.GetDataVoid());

        PyObject *pyY = PyObject_CallMethod(pyObservationManagerInstance_,
                                            function_name,
                                            format_unit, pyState);

        if (pyY == NULL)
            throw ErrorPythonUndefined("PythonObservationManager"
                                       "::ApplyTangentLinearOperator",
                                       string(function_name),
                                       "(self, state, observation)", module_);

        if (!PyArray_ISCONTIGUOUS(pyY))
            throw ErrorProcessing("PythonObservationManager"
                                  "::ApplyTangentLinearOperator",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *yArray = reinterpret_cast<PyArrayObject*>(pyY);

        int size = yArray->dimensions[0];
        y.Reallocate(size);
        memcpy(y.GetDataVoid(), yArray->data, size * sizeof(double));
    }


    //! Returns an element of the tangent linear operator.
    /*! This method is called after 'SetTime' set the time at which the
      operator is defined.
      \param[in] i row index.
      \param[in] j column index.
      \return The element (\a i, \a j) of the tangent linear operator.
    */
    double PythonObservationManager::GetTangentLinearOperator(int i, int j)
    const
    {
        char function_name[] = "GetTangentLinearOperator";
        char format_unit[] = "ii";
        PyObject *pyDouble =
            PyObject_CallMethod(pyObservationManagerInstance_,
                                function_name, format_unit, i, j);

        if (pyDouble == NULL)
            throw ErrorPythonUndefined("PythonObservationManager"
                                       "::GetTangentLinearOperator",
                                       string(function_name), "(self, i, j)",
                                       module_);

        return PyFloat_AsDouble(pyDouble);
    }


    //! Returns a row of the tangent linear operator.
    /*! This method is called after 'SetTime' set the time at which the
      operator is defined.
      \param[in] row row index.
      \param[out] tangent_operator_row the row \a row of the tangent linear
      operator.
    */
    void PythonObservationManager
    ::GetTangentLinearOperatorRow(int row,
                                  PythonObservationManager
                                  ::tangent_linear_operator_row&
                                  tangent_operator_row)
    const
    {
        char function_name[] = "GetTangentLinearOperatorRow";
        char format_unit[] = "i";

        PyObject *pyRow =
            PyObject_CallMethod(pyObservationManagerInstance_,
                                function_name, format_unit, row);

        if (pyRow == NULL)
            throw ErrorPythonUndefined("PythonObservationManager"
                                       "::GetTangentLinearOperatorRow",
                                       string(function_name),
                                       "(self, row, tangent_operator_row)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyRow))
            throw ErrorProcessing("PythonObservationManager"
                                  "::ApplyTangentLinearOperator",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *rowArray = reinterpret_cast<PyArrayObject*>(pyRow);

        int size = rowArray->dimensions[0];
        tangent_operator_row.Reallocate(size);
        memcpy(tangent_operator_row.GetDataVoid(), rowArray->data,
               size * sizeof(double));
    }


    //! Returns the tangent linear operator.
    /*! This method is called after 'SetTime' set the time at which the
      operator is defined.
      \return The matrix of the tangent linear operator.
    */
    const PythonObservationManager::tangent_linear_operator&
    PythonObservationManager::GetTangentLinearOperator()
    {
        tangent_operator_matrix_.Nullify();

        char function_name[] = "GetTangentLinearOperator";
        PyObject *pyValue = PyObject_CallMethod(pyObservationManagerInstance_,
                                                function_name, NULL);

        if (pyValue == NULL)
            throw ErrorPythonUndefined("PythonObservationManager"
                                       "::GetTangentLinearOperator",
                                       string(function_name), "(self)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyValue))
            throw ErrorProcessing("PythonModel::GetTangentLinearOperator",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *tangentArray =
            reinterpret_cast<PyArrayObject*>(pyValue);

        int size_X = tangentArray->dimensions[0];
        int size_Y = tangentArray->dimensions[1];

        tangent_operator_matrix_
            .SetData(size_X, size_Y,
                     reinterpret_cast<double*>(tangentArray->data));

        return tangent_operator_matrix_;
    }


    //! Applies the adjoint operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the operator at \a x. It is resized if
      needed.
    */
    template <class state>
    void PythonObservationManager
    ::ApplyAdjointOperator(const state& x, observation& y) const
    {
        char function_name[] = "ApplyAdjointOperator";
        char format_unit[] = "O";
        npy_intp dim[1];
        dim[0] = x.GetLength();

        PyObject *pyState = PyArray_SimpleNewFromData(1, dim, NPY_DOUBLE,
                                                      x.GetDataVoid());

        PyObject *pyY = PyObject_CallMethod(pyObservationManagerInstance_,
                                            function_name,
                                            format_unit, pyState);

        if (pyY == NULL)
            throw ErrorPythonUndefined("PythonObservationManager"
                                       "::ApplyAdjointOperator",
                                       string(function_name),
                                       "(self, state, observation)", module_);

        if (!PyArray_ISCONTIGUOUS(pyY))
            throw ErrorProcessing("PythonObservationManager"
                                  "::ApplyAdjointOperator",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *yArray = reinterpret_cast<PyArrayObject*>(pyY);

        int size = yArray->dimensions[0];
        y.Reallocate(size);
        memcpy(y.GetDataVoid(), yArray->data, size * sizeof(double));
    }


    //! Return an observation error covariance.
    /*!
      \param[in] i row index.
      \param[in] j column index.
      \return The element (\a i, \a j) of the observation error variance.
    */
    double PythonObservationManager::GetErrorVariance(int i, int j) const
    {
        char function_name[] = "GetErrorVariance";
        char format_unit[] = "ii";
        PyObject *pyDouble
            = PyObject_CallMethod(pyObservationManagerInstance_,
                                  function_name, format_unit, i, j);

        if (pyDouble == NULL)
            throw ErrorPythonUndefined("PythonObservationManager"
                                       "::GetErrorVariance",
                                       string(function_name),
                                       "(self, i, j)", module_);

        return PyFloat_AsDouble(pyDouble);
    }


    //! Returns the observation error variance.
    /*!
      \return The observation error covariance matrix.
    */
    const PythonObservationManager::error_variance&
    PythonObservationManager::GetErrorVariance()
    {
        error_variance_.Nullify();

        char function_name[] = "GetErrorVariance";
        PyObject *pyValue = PyObject_CallMethod(pyObservationManagerInstance_,
                                                function_name, NULL);

        if (pyValue == NULL)
            throw ErrorPythonUndefined("PythonModel::GetErrorVariance",
                                       string(function_name), "(self)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyValue))
            throw ErrorProcessing("PythonModel::GetErrorVariance",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *errorArray
            = reinterpret_cast<PyArrayObject*>(pyValue);

        int size_X = errorArray->dimensions[0];
        int size_Y = errorArray->dimensions[1];

        error_variance_
            .SetData(size_X, size_Y,
                     reinterpret_cast<double*>(errorArray->data));

        return error_variance_;
    }


    //! Returns the inverse of the observation error covariance matrix.
    /*!
      \return The inverse of the matrix of the observation error covariance.
    */
    const PythonObservationManager::error_variance&
    PythonObservationManager::GetErrorVarianceInverse()
    {
        error_variance_inverse_.Nullify();

        char function_name[] = "GetErrorVarianceInverse";
        PyObject *pyValue = PyObject_CallMethod(pyObservationManagerInstance_,
                                                function_name, NULL);

        if (pyValue == NULL)
            throw ErrorPythonUndefined("PythonModel"
                                       "::GetErrorVarianceInverse",
                                       string(function_name), "(self)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyValue))
            throw ErrorProcessing("PythonModel::GetErrorVarianceInverse",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *errorArray
            = reinterpret_cast<PyArrayObject*>(pyValue);

        int size_X = errorArray->dimensions[0];
        int size_Y = errorArray->dimensions[1];

        error_variance_inverse_
            .SetData(size_X, size_Y,
                     reinterpret_cast<double*>(errorArray->data));

        return error_variance_inverse_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    string PythonObservationManager::GetName() const
    {
        char function_name[] = "GetName";
        PyObject *pyName = PyObject_CallMethod(pyObservationManagerInstance_,
                                               function_name, NULL);
        if (pyName == NULL)
            ErrorPythonUndefined("PythonObseravationManager::GetName",
                                 string(function_name), "(self)", module_);

        return PyString_AsString(pyName);
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    void PythonObservationManager
    ::Message(string message)
    {
        char function_name[] = "Message";
        char format_unit[] = "s";

        if (is_module_initialized_ &&
            PyObject_CallMethod(pyObservationManagerInstance_, function_name,
                                format_unit, message.c_str()) != Py_None)
            throw ErrorPythonUndefined("PythonObservationManager::Message",
                                       string(function_name),
                                       "(self, string)", module_);
    }


    ////////////////////
    // PRIVATE METHOD //
    ////////////////////


    //! Returns the error message for non-contiguous Python arrays.
    /*
      \param[in] function_name name of the function that returned the
      non-contiguous array.
      \return The error message for non-contiguous Python arrays.
    */
    string PythonObservationManager
    ::ErrorMessageNotContiguous(string function_name) const
    {
        return "An output array of " + function_name
            + " is not contiguous in memory or not in C-style order.";
    }


} // namespace Verdandi

#define VERDANDI_FILE_OBSERVATION_MANAGER_PYTHONOBSERVATIONMANAGER_CXX
#endif
