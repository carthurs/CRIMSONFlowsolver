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


#ifndef VERDANDI_FILE_MODEL_PYTHONMODEL_CXX


#include "PythonModel.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Constructor.
    PythonModel::PythonModel()
    {
        Py_Initialize();
        import_array();
        is_module_initialized_ = false;
    }


    //! Constructor.
    /*! It reads the initial condition and the time settings.
      \param[in] configuration_file path to the configuration file.
    */
    PythonModel::PythonModel(string configuration_file)
    {
        Py_Initialize();
        import_array();
        is_module_initialized_ = false;

        Initialize(configuration_file);
    }


    //! Destructor.
    PythonModel::~PythonModel()
    {
        Py_Finalize();

        state_.Nullify();
        uncertain_parameter_.Nullify();
        parameter_correlation_.Nullify();
        parameter_variance_.Nullify();
        parameter_parameter_.Nullify();
        state_error_variance_.Nullify();
        state_error_variance_inverse_.Nullify();
    }


    /////////////////////
    // INITIALIZATIONS //
    /////////////////////


    //! Initializes the model.
    /*! It reads the initial condition and the time settings.
      \param[in] configuration_file configuration file.
    */
    void PythonModel::Initialize(string configuration_file)
    {
        VerdandiOps configuration(configuration_file);

        configuration.Set("python_model.module", module_);
        configuration.Set("python_model.directory", directory_);
        configuration.Set("python_model.class_name", class_name_);

        pyModelFile_ = PyString_FromString(module_.c_str());

        if (directory_ != "")
        {
            string PythonPath = "sys.path.append(\"" + directory_ + "\")";
            PyRun_SimpleString("import sys");
            PyRun_SimpleString(PythonPath.c_str());
        }

        pyModelModule_ = PyImport_Import(pyModelFile_);

        pyModelDict_ = PyModule_GetDict(pyModelModule_);

        pyModelClass_ = PyDict_GetItemString(pyModelDict_,
                                             class_name_.c_str());

        PyObject *file_path =
            PyString_FromString(configuration.GetFilePath().c_str());
        PyObject *args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, file_path);

        if (!PyCallable_Check(pyModelClass_))
            throw ErrorConfiguration("PythonModel::PythonModel",
                                     "The class \"" + class_name_
                                     + "\" is not defined in the module \""
                                     + module_ + "\".");

        pyModelInstance_ = PyObject_CallObject(pyModelClass_, args);

        is_module_initialized_ = true;

        Py_DECREF(pyModelFile_);
        Py_DECREF(pyModelModule_);
    }


    //! Initializes the current time step for the model.
    void PythonModel::InitializeStep()
    {
        char function_name[] = "InitializeStep";
        if (PyObject_CallMethod(pyModelInstance_, function_name, NULL)
            != Py_None)
            throw ErrorPythonUndefined("PythonModel::InitializeStep",
                                       string(function_name), "(self)",
                                       module_);
    }


    ////////////////
    // PROCESSING //
    ////////////////


    //! Advances one step forward in time.
    /*! \f[x^f_{h+1} = \mathcal{M}_h(x^a_h, p_h)\,.\f] */
    void PythonModel::Forward()
    {
        char function_name[] = "Forward";
        if (PyObject_CallMethod(pyModelInstance_, function_name, NULL)
            != Py_None)
            throw ErrorPythonUndefined("PythonModel::Forward",
                                       string(function_name), "(self)",
                                       module_);
    }


    //! Performs one step backward in adjoint model.
    /*!
      \param[in] observation_term \f$ H^T R^{-1}(y - Hx) \f$.
    */
    void PythonModel::BackwardAdjoint(state& observation_term)
    {
        char function_name[] = "BackwardAdjoint";
        char format_unit[] = "O";
        npy_intp dim[1];
        dim[0] = observation_term.GetLength();

        PyObject *pyState =
            PyArray_SimpleNewFromData(1, dim, NPY_DOUBLE,
                                      observation_term.GetDataVoid());

        if (PyObject_CallMethod(pyModelInstance_, function_name, format_unit,
                                pyState) != Py_None)
            throw ErrorPythonUndefined("PythonModel::BackwardAdjoint",
                                       string(function_name), "(self)",
                                       module_);
    }


    //! Checks whether the model has finished.
    /*!
      \return True if the simulation is done, false otherwise.
    */
    bool PythonModel::HasFinished() const
    {
        char function_name[] = "HasFinished";
        PyObject *has_finished = PyObject_CallMethod(pyModelInstance_,
                                                     function_name, NULL);
        if (has_finished == NULL)
            throw ErrorPythonUndefined("PythonModel::HasFinished",
                                       string(function_name), "(self)",
                                       module_);

        return PyInt_AsLong(has_finished) != 0;
    }


    //! Saves the simulated data.
    void PythonModel::Save()
    {
    }


    //! Finalizes the current time step for the model.
    void PythonModel::FinalizeStep()
    {
        char function_name[] = "FinalizeStep";
        if (PyObject_CallMethod(pyModelInstance_, function_name, NULL)
            != Py_None)
            throw ErrorPythonUndefined("PythonModel::FinalizeStep",
                                       string(function_name), "(self)",
                                       module_);
    }


    //! Finalizes the model.
    void PythonModel::Finalize()
    {
        char function_name[] = "Finalize";
        if (PyObject_CallMethod(pyModelInstance_, function_name, NULL)
            != Py_None)
            throw ErrorPythonUndefined("PythonModel::Finalize",
                                       string(function_name), "(self)",
                                       module_);
    }


    ///////////////
    // OPERATORS //
    ///////////////


    //! Applies the model to a given vector.
    /*! The current state of the model is modified.
      \param[in,out] x a vector.
      \param[in] forward Boolean to indicate if the model has to go on to the
      next step.
      \param[in] preserve_state Boolean to indicate if the model state has to
      be preserved.
    */
    void PythonModel::ApplyOperator(state& x,
                                    bool forward, bool preserve_state)
    {
        char function_name[] = "ApplyOperator";
        char format_unit[] = "Oii";
        npy_intp dim[1];
        dim[0] = x.GetLength();

        PyObject *pyState =  PyArray_SimpleNewFromData(1, dim, NPY_DOUBLE,
                                                       x.GetDataVoid());
        if (PyObject_CallMethod(pyModelInstance_, function_name, format_unit,
                                pyState, forward, preserve_state) != Py_None)
            throw ErrorPythonUndefined("PythonModel::ApplyOperator",
                                       string(function_name),
                                       "(self, state, forward,"
                                       " preserve_state)",
                                       module_);
    }


    //! Applies the tangent linear model to a given vector.
    /*!
      \param[in,out] x on entry, a vector to which the tangent linear model
      should be applied; on exit, the result.
    */
    void PythonModel::ApplyTangentLinearOperator(state& x)
    {
        char function_name[] = "ApplyTangentLinearOperator";
        char format_unit[] = "O";
        npy_intp dim[1];
        dim[0] = x.GetLength();

        PyObject *pyState =  PyArray_SimpleNewFromData(1, dim, NPY_DOUBLE,
                                                       x.GetDataVoid());

        if (PyObject_CallMethod(pyModelInstance_, function_name, format_unit,
                                pyState) != Py_None)
            throw ErrorPythonUndefined("PythonModel"
                                       "::ApplyTangentLinearOperator",
                                       string(function_name), "(self, state)",
                                       module_);
    }


    //! Gets the tangent linear model.
    /*!
      \param[out] A the matrix of the tangent linear model.
    */
    void PythonModel::
    GetTangentLinearOperator(tangent_linear_operator& A) const
    {
        char function_name[] = "GetTangentLinearOperator";

        PyObject *pyValue = PyObject_CallMethod(pyModelInstance_,
                                                function_name, NULL);
        if (pyValue == NULL)
            throw ErrorPythonUndefined("PythonModel"
                                       "::GetTangentLinearOperator",
                                       string(function_name), "(self)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyValue))
            throw ErrorProcessing("PythonModel::GetTangentLinearOperator",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *stateArray =
            reinterpret_cast<PyArrayObject*>(pyValue);

        int size_X = (stateArray->dimensions)[0];
        int size_Y = (stateArray->dimensions)[1];
        A.Reallocate(size_X, size_Y);
        memcpy(A.GetDataVoid(), stateArray->data,
               size_X * size_Y * sizeof(double));
    }


    ////////////////////
    // ACCESS METHODS //
    ////////////////////


    //! Returns the current time.
    /*!
      \return The current time.
    */
    double PythonModel::GetTime() const
    {
        char function_name[] = "GetTime";
        PyObject *pyTime = PyObject_CallMethod(pyModelInstance_,
                                               function_name, NULL);
        if (pyTime == NULL)
            throw ErrorPythonUndefined("PythonModel::GetTime",
                                       string(function_name), "(self)",
                                       module_);

        return PyFloat_AsDouble(pyTime);
    }


    //! Sets the current time.
    /*!
      \return The time step.
    */
    void PythonModel::SetTime(double time)
    {
        char function_name[] = "SetTime";
        char format_unit[] = "d";
        if (PyObject_CallMethod(pyModelInstance_, function_name,
                                format_unit, time) != Py_None)
            throw ErrorPythonUndefined("PythonModel::SetTime",
                                       string(function_name), "(self, time)",
                                       module_);
    }


    //! Returns the state vector size.
    /*!
      \return The state vector size.
    */
    int PythonModel::GetNstate() const
    {
        char function_name[] = "GetNstate";
        PyObject *Nstate = PyObject_CallMethod(pyModelInstance_,
                                               function_name, NULL);
        if (Nstate == NULL)
            throw ErrorPythonUndefined("PythonModel::GetNstate",
                                       string(function_name), "(self)",
                                       module_);

        return PyInt_AsLong(Nstate);
    }


    //! Returns the full state vector size.
    /*!
      \return The full state vector size.
    */
    int PythonModel::GetNfull_state() const
    {
        char function_name[] = "GetNFullstate";
        PyObject *Nfull_state = PyObject_CallMethod(pyModelInstance_,
                                                    function_name, NULL);
        if (Nfull_state == NULL)
            throw ErrorPythonUndefined("PythonModel::GetNfull_state",
                                       string(function_name), "(self)",
                                       module_);

        return PyInt_AsLong(Nfull_state);
    }


    //! Provides the controlled state vector.
    /*!
      \return the controlled state vector.
    */
    PythonModel::state& PythonModel::GetState()
    {
        state_.Nullify();

        char function_name[] = "GetState";

        PyObject *pyValue = PyObject_CallMethod(pyModelInstance_,
                                                function_name, NULL);

        if (pyValue == NULL)
            throw ErrorPythonUndefined("PythonModel::GetState",
                                       string(function_name), "(self)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyValue))
            throw ErrorProcessing("PythonModel::GetState",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *stateArray = reinterpret_cast<PyArrayObject*>(pyValue);

        int size = (stateArray->dimensions)[0];

        state_.SetData(size,
                       reinterpret_cast<double*> (stateArray->data));

        return state_;
    }


    //! Performs some calculations when the update of the model state is done.
    void PythonModel::StateUpdated()
    {
        char function_name[] = "StateUpdated";

        if (PyObject_CallMethod(pyModelInstance_, function_name, NULL)
            != Py_None)
            throw ErrorPythonUndefined("PythonModel::StateUpdated",
                                       string(function_name), "(self)",
                                       module_);
    }


    //! Provides the state lower bound.
    /*!
      \return the state lower bound (componentwise).
    */
    PythonModel::state& PythonModel::GetStateLowerBound()
    {
        char function_name[] = "GetStateLowerBound";

        PyObject *pyValue = PyObject_CallMethod(pyModelInstance_,
                                                function_name, NULL);

        if (pyValue == NULL)
            throw ErrorPythonUndefined("PythonModel::GetStateLowerBound",
                                       string(function_name), "(self)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyValue))
            throw ErrorProcessing("PythonModel::GetStateLowerBound",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *lower_bound_array =
            reinterpret_cast<PyArrayObject*>(pyValue);

        int size = (lower_bound_array->dimensions)[0];

        lower_bound_.SetData(size,reinterpret_cast<double*>
                             (lower_bound_array->data));

        return lower_bound_;
    }


    //! Provides the state upper bound.
    /*!
      \return the state upper bound (componentwise).
    */
    PythonModel::state& PythonModel::GetStateUpperBound()
    {
        char function_name[] = "GetStateUpperBound";

        PyObject *pyValue = PyObject_CallMethod(pyModelInstance_,
                                                function_name, NULL);

        if (pyValue == NULL)
            throw ErrorPythonUndefined("PythonModel::GetStateUpperBound",
                                       string(function_name), "(self)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyValue))
            throw ErrorProcessing("PythonModel::GetStateUpperBound",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *upper_bound_array =
            reinterpret_cast<PyArrayObject*>(pyValue);

        int size = (upper_bound_array->dimensions)[0];
        upper_bound_.SetData(size, reinterpret_cast<double*>
                             (upper_bound_array->data));

        return upper_bound_;
    }


    //! Provides the full state vector.
    /*!
      \param[out] x the full state vector.
    */
    PythonModel::state& PythonModel::GetFullState()
    {
        full_state_.Nullify();

        char function_name[] = "GetState";

        PyObject *pyValue = PyObject_CallMethod(pyModelInstance_,
                                                function_name, NULL);

        if (pyValue == NULL)
            throw ErrorPythonUndefined("PythonModel::GetFullState",
                                       string(function_name), "(self)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyValue))
            throw ErrorProcessing("PythonModel::GetState",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *stateArray = reinterpret_cast<PyArrayObject*>(pyValue);

        int size = (stateArray->dimensions)[0];

        full_state_.SetData(size,
                            reinterpret_cast<double*> (stateArray->data));

        return full_state_;
    }


    //! Performs some calculations when the update of the model state is done.
    void PythonModel::FullStateUpdated()
    {
        char function_name[] = "FullStateUpdated";

        if (PyObject_CallMethod(pyModelInstance_, function_name, NULL)
            != Py_None)
            throw ErrorPythonUndefined("PythonModel::FullStateUpdated",
                                       string(function_name), "(self)",
                                       module_);
    }


    //! Returns the adjoint state vector.
    /*!
      \param[out] state_adjoint the adjoint state vector.
    */
    void PythonModel::GetAdjointState(state& state_adjoint)
    {
        char function_name[] = "GetAdjointState";

        PyObject *pyValue = PyObject_CallMethod(pyModelInstance_,
                                                function_name, NULL);

        if (pyValue == NULL)
            throw ErrorPythonUndefined("PythonModel::GetAdjointState",
                                       string(function_name), "(self)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyValue))
            throw ErrorProcessing("PythonModel::GetAdjointState",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *stateArray = reinterpret_cast<PyArrayObject*>(pyValue);

        int size = (stateArray->dimensions)[0];
        state_adjoint.Reallocate(size);
        memcpy(state_adjoint.GetDataVoid(), stateArray->data,
               size*sizeof(double));
    }


    //! Sets the adjoint state vector.
    /*!
      \param[out] state_adjoint the adjoint state vector.
    */
    void PythonModel::SetAdjointState(const state& state_adjoint)
    {
        char function_name[] = "SetAdjointState";
        char format_unit[] = "O";
        npy_intp dim[1];
        dim[0] = state_adjoint.GetLength();

        PyObject *pyState =
            PyArray_SimpleNewFromData(1, dim, NPY_DOUBLE,
                                      state_adjoint.GetDataVoid());

        if (PyObject_CallMethod(pyModelInstance_, function_name,
                                format_unit, pyState) != Py_None)
            throw ErrorPythonUndefined("PythonModel::SetAdjointState",
                                       string(function_name),
                                       "(self, state_adjoint)", module_);
    }


    //! Returns the number of parameters to be perturbed.
    /*!
      \return The number of parameters to be perturbed.
    */
    int PythonModel::GetNparameter()
    {
        char function_name[] = "GetNparameter";
        PyObject *Nparameter = PyObject_CallMethod(pyModelInstance_,
                                                   function_name, NULL);

        if (Nparameter == NULL)
            throw ErrorPythonUndefined("PythonModel::GetNparameter",
                                       string(function_name), "(self)",
                                       module_);

        return PyInt_AsLong(Nparameter);
    }


    //! Returns the i-th uncertain parameter.
    /*!
      \param[in] i index of the parameter.
      \return The vector associated with the i-th parameter.
    */
    PythonModel::uncertain_parameter& PythonModel::GetParameter(int i)
    {
        uncertain_parameter_.Nullify();

        char function_name[] = "GetParameter";
        char format_unit[] = "i";
        PyObject *pyParameter = PyObject_CallMethod(pyModelInstance_,
                                                    function_name,
                                                    format_unit, i);

        if (pyParameter == NULL)
            throw ErrorPythonUndefined("PythonModel::GetParameter",
                                       string(function_name), "(self, i)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyParameter))
            throw ErrorProcessing("PythonModel::GetParameter",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *parameterArray =
            reinterpret_cast<PyArrayObject*>(pyParameter);

        int size = (parameterArray->dimensions)[0];

        uncertain_parameter_.SetData(size,
                                     reinterpret_cast<double*>
                                     (parameterArray->data));

        return uncertain_parameter_;
    }


    //! Sets the i-th uncertain parameter.
    /*!
      \param[in] i index of the parameter.
      \param[in] parameter The new uncertain parameter.
    */
    void PythonModel::SetParameter(int i, uncertain_parameter& parameter)
    {
        char function_name[] = "SetParameter";
        char format_unit[] = "iO";
        npy_intp dim[1];
        dim[0] = parameter.GetLength();

        PyObject *pyParameter =
            PyArray_SimpleNewFromData(1, dim, NPY_DOUBLE,
                                      parameter.GetDataVoid());

        if (PyObject_CallMethod(pyModelInstance_, function_name,
                                format_unit, i, pyParameter) != Py_None)
            throw ErrorPythonUndefined("PythonModel::SetParameter",
                                       string(function_name),
                                       "(self, i, parameter)", module_);
    }


    //! Returns the correlation between the uncertain parameters.
    /*!
      \param[in] i index of the parameter.
      \return The correlation between the uncertain parameters.
    */
    Vector<double>& PythonModel::GetParameterCorrelation(int i)
    {
        parameter_correlation_.Nullify();

        char function_name[] = "GetParameterCorrelation";
        char format_unit[] = "i";
        PyObject *pyCorrelation = PyObject_CallMethod(pyModelInstance_,
                                                      function_name,
                                                      format_unit, i);

        if (pyCorrelation == Py_None) return parameter_correlation_;

        if (pyCorrelation == NULL)
            throw ErrorPythonUndefined("PythonModel::GetParameterCorrelation",
                                       string(function_name), "(self, i)",
                                       module_);


        if (!PyArray_ISCONTIGUOUS(pyCorrelation))
            throw ErrorProcessing("PythonModel::GetParameterCorrelation",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *correlationArray =
            reinterpret_cast<PyArrayObject*>(pyCorrelation);

        int size = (correlationArray->dimensions)[0];

        parameter_correlation_.SetData(size,
                                       reinterpret_cast<double*>
                                       (correlationArray->data));

        return parameter_correlation_;
    }


    //! Returns the PDF of the i-th parameter.
    /*!
      \param[in] i uncertain-variable index.
      \return The PDF of the i-th parameter.
    */
    string PythonModel::GetParameterPDF(int i)
    {
        char function_name[] = "GetParameterPDF";
        char format_unit[] = "i";
        PyObject *pyPDF = PyObject_CallMethod(pyModelInstance_,
                                              function_name, format_unit, i);
        if (pyPDF == NULL)
            throw ErrorPythonUndefined("PythonModel::GetParameterPDF",
                                       string(function_name), "(self, i)",
                                       module_);

        return PyString_AsString(pyPDF);
    }


    /*! \brief Returns the covariance matrix associated with the i-th
      parameter.*/
    /*!
      \param[in] i parameter index.
      \return The covariance matrix associated with the i-th parameter.
    */
    Matrix<double>& PythonModel::GetParameterVariance(int i)
    {
        parameter_variance_.Nullify();

        char function_name[] = "GetParameterVariance";
        char format_unit[] = "i";
        PyObject *pyVariance = PyObject_CallMethod(pyModelInstance_,
                                                   function_name,
                                                   format_unit, i);

        if (pyVariance == NULL)
            throw ErrorPythonUndefined("PythonModel::GetParameterVariance",
                                       string(function_name), "(self, i)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyVariance))
            throw ErrorProcessing("PythonModel::GetParameterVariance",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *varianceArray =
            reinterpret_cast<PyArrayObject*>(pyVariance);

        int size_X = (varianceArray->dimensions)[0];
        int size_Y = (varianceArray->dimensions)[1];

        parameter_variance_.SetData(size_X, size_Y,
                                    reinterpret_cast<double*>
                                    (varianceArray->data));

        return parameter_variance_;
    }


    //! Returns parameters associated with the PDF of some model parameter.
    /*! In case of normal or log-normal distribution, the parameters are
      clipping parameters.
      \param[in] i model parameter index.
      \return The parameters associated with the i-th parameter.
    */
    Vector<double>& PythonModel::GetParameterParameter(int i)
    {
        parameter_parameter_.Nullify();

        char function_name[] = "GetParameterParameter";
        char format_unit[] = "i";
        PyObject *pyParameter = PyObject_CallMethod(pyModelInstance_,
                                                    function_name,
                                                    format_unit, i);

        if (pyParameter == NULL)
            throw ErrorPythonUndefined("PythonModel::GetParameterParameter",
                                       string(function_name), "(self, i)",
                                       module_);

        if (pyParameter == Py_None) return parameter_parameter_;

        if (!PyArray_ISCONTIGUOUS(pyParameter))
            throw ErrorProcessing("PythonModel::GetParameterParameter",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *parameterArray =
            reinterpret_cast<PyArrayObject*>(pyParameter);

        int size = (parameterArray->dimensions)[0];

        parameter_parameter_.SetData(size,
                                     reinterpret_cast<double*>
                                     (parameterArray->data));

        return parameter_parameter_;
    }


    //! Returns the perturbation option of the i-th parameter.
    /*!
      \param[in] i parameter index.
      \return The perturbation option of the i-th parameter.
    */
    string PythonModel::GetParameterOption(int i)
    {
        char function_name[] = "GetParameterOption";
        char format_unit[] = "i";
        PyObject *pyOption = PyObject_CallMethod(pyModelInstance_,
                                                 function_name,
                                                 format_unit, i);
        if (pyOption == NULL)
            throw ErrorPythonUndefined("PythonModel::GetParameterOption",
                                       string(function_name), "(self, i)",
                                       module_);

        return PyString_AsString(pyOption);
    }


    ////////////
    // ERRORS //
    ////////////


    //! Computes a row of the variance of the state error.
    /*!
      \param[in] row row index.
      \param[out] P_row the row with index \a row in the state error variance.
    */
    void PythonModel
    ::GetStateErrorVarianceRow(int row, state_error_variance_row& P_row)
    {
        char function_name[] = "GetStateErrorVarianceRow";
        char format_unit[] = "i";

        PyObject *pyValue = PyObject_CallMethod(pyModelInstance_,
                                                function_name, format_unit,
                                                row);

        if (pyValue == NULL)
            throw ErrorPythonUndefined("PythonModel"
                                       "::GetStateErrorVarianceRow",
                                       string(function_name), "(self, row)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyValue))
            throw ErrorProcessing("PythonModel::GetStateErrorVarianceRow",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *stateArray = reinterpret_cast<PyArrayObject*>(pyValue);

        int size = (stateArray->dimensions)[0];
        P_row.Reallocate(size);
        memcpy(P_row.GetDataVoid(), stateArray->data,
               size*sizeof(double));
    }


    //! Returns the state error variance.
    /*!
      \return The state error variance.
    */
    PythonModel::state_error_variance&
    PythonModel::GetStateErrorVariance()
    {
        state_error_variance_.Nullify();

        char function_name[] = "GetStateErrorVariance";
        PyObject *pyValue = PyObject_CallMethod(pyModelInstance_,
                                                function_name, NULL);

        if (pyValue == NULL)
            throw ErrorPythonUndefined("PythonModel::GetStateErrorVariance",
                                       string(function_name), "(self)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyValue))
            throw ErrorProcessing("PythonModel::GetStateErrorVariance",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *stateArray = reinterpret_cast<PyArrayObject*>(pyValue);

        int size_X = (stateArray->dimensions)[0];
        int size_Y = (stateArray->dimensions)[1];

        state_error_variance_.SetData(size_X, size_Y,
                                      reinterpret_cast<double*>
                                      (stateArray->data));

        return state_error_variance_;
    }


    /*! Returns a decomposition of the state error covariance matrix (\f$B\f$)
      as a product \f$LUL^T\f$.
    */
    /*!
      \param[out] L the matrix \f$L\f$.
      \param[out] U the matrix \f$U\f$.
    */
    void PythonModel::GetStateErrorVarianceSqrt(state_error_variance& L,
                                                state_error_variance& U)
    {

        char function_name[] = "GetStateErrorVarianceSqrt";
        PyObject *pyValue = PyObject_CallMethod(pyModelInstance_,
                                                function_name, NULL);

        if (pyValue == NULL)
            throw ErrorPythonUndefined("PythonModel"
                                       "::GetStateErrorVarianceSqrt",
                                       string(function_name), "(self)",
                                       module_);

        PyObject *py_L = PyTuple_GetItem(pyValue, (Py_ssize_t) 0);
        PyObject *py_U = PyTuple_GetItem(pyValue, (Py_ssize_t) 1);

        if (!PyArray_ISCONTIGUOUS(py_L))
            throw ErrorProcessing("PythonModel::GetStateErrorVarianceSqrt",
                                  ErrorMessageNotContiguous(function_name));
        if (!PyArray_ISCONTIGUOUS(py_U))
            throw ErrorProcessing("PythonModel::GetStateErrorVarianceSqrt",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *L_array = reinterpret_cast<PyArrayObject*>(py_L);

        int size_x = (L_array->dimensions)[0];
        int size_y = (L_array->dimensions)[1];
        L.Reallocate(size_x, size_y);
        memcpy(L.GetDataVoid(), L_array->data,
               size_x*size_y*sizeof(double));

        PyArrayObject *U_array = reinterpret_cast<PyArrayObject*>(py_U);

        size_x = (U_array->dimensions)[0];
        size_y = (U_array->dimensions)[1];
        U.Reallocate(size_x, size_y);
        memcpy(U.GetDataVoid(), U_array->data,
               size_x*size_y*sizeof(double));
    }


    //! Returns the inverse of the background error variance (\f$B^{-1}\f$).
    /*!
      \return The inverse of the background error variance (\f$B^{-1}\f$).
    */
    const PythonModel::state_error_variance&
    PythonModel::GetStateErrorVarianceInverse()
    {
        state_error_variance_inverse_.Nullify();

        char function_name[] = "GetStateErrorVarianceInverse";
        PyObject *pyValue = PyObject_CallMethod(pyModelInstance_,
                                                function_name, NULL);

        if (pyValue == NULL)
            throw ErrorPythonUndefined("PythonModel"
                                       "::GetStateErrorVarianceInverse",
                                       string(function_name), "(self)",
                                       module_);

        if (!PyArray_ISCONTIGUOUS(pyValue))
            throw ErrorProcessing("PythonModel::GetStateErrorVarianceInverse",
                                  ErrorMessageNotContiguous(function_name));

        PyArrayObject *stateArray = reinterpret_cast<PyArrayObject*>(pyValue);

        int size_X = (stateArray->dimensions)[0];
        int size_Y = (stateArray->dimensions)[1];

        state_error_variance_inverse_.SetData(size_X, size_Y,
                                              reinterpret_cast<double*>
                                              (stateArray->data));

        return state_error_variance_inverse_;
    }


    //! Returns the name of the model.
    /*!
      \return The name of the model.
    */
    string PythonModel::GetName() const
    {
        char function_name[] = "GetName";
        PyObject *pyName = PyObject_CallMethod(pyModelInstance_,
                                               function_name, NULL);
        if (pyName == NULL)
            ErrorPythonUndefined("PythonModel::GetName",
                                 string(function_name), "(self)", module_);

        return PyString_AsString(pyName);
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    void PythonModel::Message(string message)
    {
        char function_name[] = "Message";
        char format_unit[] = "s";

        if (is_module_initialized_ &&
            PyObject_CallMethod(pyModelInstance_, function_name,
                                format_unit, message.c_str()) != Py_None)
            throw ErrorPythonUndefined("PythonModel::Message",
                                       string(function_name),
                                       "(self, string)", module_);
    }


    /////////////////////
    // PRIVATE METHODS //
    /////////////////////


    //! Returns the error message for non-contiguous Python arrays.
    /*
      \param[in] function_name name of the function that returned the
      non-contiguous array.
      \return The error message for non-contiguous Python arrays.
    */
    string PythonModel::ErrorMessageNotContiguous(string function_name) const
    {
        return "An output array of " + function_name
            + " is not contiguous in memory or not in C-style order.";
    }


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_PYTHONMODEL_CXX
#endif
