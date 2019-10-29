// Copyright (C) 2010 INRIA
// Author(s): Marc Fragu, Vivien Mallet, Claire Mouton
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


#ifndef VERDANDI_FILE_OBSERVATIONMANAGER_OBSERVATIONMANAGERTEMPLATE_CXX


#include "ObservationManagerTemplate.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    ObservationManagerTemplate::ObservationManagerTemplate()
    {
        throw ErrorUndefined(
            "ObservationManagerTemplate::ObservationManagerTemplate()");
    }


    //! Destructor.
    ObservationManagerTemplate::~ObservationManagerTemplate()
    {
        // Operations to be performed when the object is destroyed.
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
    void ObservationManagerTemplate
    ::Initialize(const Model& model, string configuration_file)
    {
        throw ErrorUndefined(
            "void ObservationManagerTemplate"
            "::Initialize(const Model& model, string configuration_file)");
    }


    //! Activates or deactivates the option 'discard_observation'.
    /*!
      \param[in] discard_observation if set to true, each observation will be
      used at most one time.
    */
    void ObservationManagerTemplate
    ::DiscardObservation(bool discard_observation)
    {
        throw ErrorUndefined(
            "void ObservationManagerTemplate"
            "::DiscardObservation(bool discard_observation)");
    }


    //! Sets the time of observations to be loaded.
    /*!
      \param[in] model the model.
      \param[in] time a given time.
    */
    template <class Model>
    void ObservationManagerTemplate
    ::SetTime(const Model& model, double time)
    {
        throw ErrorUndefined(
            "void ObservationManagerTemplate"
            "::SetTime(const Model& model, double time)");
    }


    /////////////////
    // OBSERVATION //
    /////////////////


    //! Returns the observations.
    /*! This method is called after 'SetTime' set the time at which the
      observations are requested.
      \param[out] observation observation vector.
    */
    void ObservationManagerTemplate::GetObservation(
        ObservationManagerTemplate::observation& observation)
    {
        throw ErrorUndefined(
            "void ObservationManagerTemplate::GetObservation("
            "ObservationManagerTemplate::observation& observation)");
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
    void ObservationManagerTemplate
    ::GetInnovation(const state& x,
                    ObservationManagerTemplate::observation& innovation)
    {
        throw ErrorUndefined(
            "void ObservationManagerTemplate::GetInnovation(const state& x,"
            "ObservationManagerTemplate::observation& innovation)");
    }


    ////////////
    // ACCESS //
    ////////////


    //! Indicates if some observations are available at a given time.
    /*!
      \param[in] time a given time.
    */
    bool ObservationManagerTemplate::HasObservation(double time)
    {
        throw ErrorUndefined(
            "bool ObservationManagerTemplate::HasObservation(double time)");
    }


    //! Indicates if some observations are available at current time.
    bool ObservationManagerTemplate::HasObservation() const
    {
        throw ErrorUndefined(
            "bool ObservationManagerTemplate::HasObservation() const");
    }


    //! Returns the number of available observations.
    /*!
      \return The total number of observation at current time.
    */
    int ObservationManagerTemplate::GetNobservation() const
    {
        throw ErrorUndefined(
            "int ObservationManagerTemplate::GetNobservation() const");
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
    void ObservationManagerTemplate
    ::ApplyOperator(const state& x, observation& y) const
    {
        throw ErrorUndefined(
            "void ObservationManagerTemplate"
            "::ApplyOperator(const state& x, observation& y) const");
    }


    //! Applies the tangent linear operator to a given vector.
    /*! This method is called after 'SetTime' set the time at which the
      operator is defined.
      \param[in] x a vector.
      \param[out] y the value of the tangent linear operator applied to \a
      x. It is resized if needed.
    */
    template <class state>
    void ObservationManagerTemplate
    ::ApplyTangentLinearOperator(const state& x, observation& y) const
    {
        throw ErrorUndefined(
            "void ObservationManagerTemplate::ApplyTangentLinearOperator"
            "(const state& x, observation& y) const");
    }


    //! Returns an element of the tangent linear operator.
    /*! This method is called after 'SetTime' set the time at which the
      operator is defined.
      \param[in] i row index.
      \param[in] j column index.
      \return The element (\a i, \a j) of the tangent linear operator.
    */
    double ObservationManagerTemplate::GetTangentLinearOperator(int i, int j)
        const
    {
        throw ErrorUndefined(
            "double ObservationManagerTemplate::"
            "GetTangentLinearOperator(int i, int j) const");
    }


    //! Returns a row of the tangent linear operator.
    /*! This method is called after 'SetTime' set the time at which the
      operator is defined.
      \param[in] row row index.
      \param[out] tangent_operator_row the row \a row of the tangent linear
      operator.
    */
    void ObservationManagerTemplate
    ::GetTangentLinearOperatorRow(int row,
                                  ObservationManagerTemplate
                                  ::tangent_linear_operator_row&
                                  tangent_operator_row)
        const
    {
        throw ErrorUndefined(
            " void ObservationManagerTemplate"
            "::GetTangentLinearOperatorRow(int row,"
            "ObservationManagerTemplate::tangent_linear_operator_row&"
            "tangent_operator_row) const");
    }


    //! Returns the tangent linear operator.
    /*! This method is called after 'SetTime' set the time at which the
      operator is defined.
      \return The matrix of the tangent linear operator.
    */
    const ObservationManagerTemplate::tangent_linear_operator&
    ObservationManagerTemplate::GetTangentLinearOperator() const
    {
        throw ErrorUndefined(
            "const typename"
            "ObservationManagerTemplate::tangent_linear_operator&"
            "ObservationManagerTemplate::GetTangentLinearOperator() const");
    }


    //! Applies the adjoint operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the operator at \a x. It is resized if
      needed.
    */
    template <class state>
    void ObservationManagerTemplate
    ::ApplyAdjointOperator(const state& x, observation& y) const
    {
        throw ErrorUndefined(
            "void ObservationManagerTemplate"
            "::ApplyAdjointOperator(const state& x, observation& y) const");
    }


    //! Return an observation error covariance.
    /*!
      \param[in] i row index.
      \param[in] j column index.
      \return The element (\a i, \a j) of the observation error variance.
    */
    double ObservationManagerTemplate::GetErrorVariance(int i, int j) const
    {
        throw ErrorUndefined(
            "double ObservationManagerTemplate"
            "::GetErrorVariance(int i, int j) const");
    }


    //! Returns the observation error variance.
    /*!
      \return The observation error covariance matrix.
    */
    const ObservationManagerTemplate::error_variance&
    ObservationManagerTemplate::GetErrorVariance() const
    {
        throw ErrorUndefined(
            "const typename"
            "ObservationManagerTemplate::error_variance& "
            "ObservationManagerTemplate::GetErrorVariance() const");
    }


    //! Returns the inverse of the observation error covariance matrix.
    /*!
      \return The inverse of the matrix of the observation error covariance.
    */
    const ObservationManagerTemplate::error_variance&
    ObservationManagerTemplate::GetErrorVarianceInverse() const
    {
        throw ErrorUndefined(
            "const typename"
            "ObservationManagerTemplate::error_variance& "
            "ObservationManagerTemplate::GetErrorVarianceInverse() const");
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    string ObservationManagerTemplate::GetName() const
    {
        return "ObservationManagerTemplate";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    void ObservationManagerTemplate
    ::Message(string message)
    {
        // Put here any processing you need.
    }


}

#define VERDANDI_FILE_OBSERVATIONMANAGER_OBSERVATIONMANAGERTEMPLATE_CXX
#endif
