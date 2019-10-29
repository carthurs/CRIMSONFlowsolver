// Copyright (C) 2008-2012 INRIA
// Author(s): Vivien Mallet
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


#ifndef VERDANDI_FILE_SEQUENTIAL_AGGREGATION_DISCOUNTEDRIDGEREGRESSION_CXX
#define VERDANDI_FILE_SEQUENTIAL_AGGREGATION_DISCOUNTEDRIDGEREGRESSION_CXX


#include "DiscountedRidgeRegression.hxx"


namespace Verdandi
{


    //! Default constructor.
    /*! The penalization is set to 0, the parameter gamma is set to 0 and the
      parameter p is set to 2. */
    template <class T>
    DiscountedRidgeRegression<T>::DiscountedRidgeRegression():
        BaseForecaster<T>(), penalization_(0), gamma_(0), p_(2)
    {
        this->parameter_.Reallocate(3);
        this->parameter_(0) = T(0);
        this->parameter_(1) = T(0);
        this->parameter_(2) = T(2);
    }


    //! Main constructor.
    /*!
      \param[in] penalization penalization.
      \param[in] gamma parameter gamma.
      \param[in] p parameter p.
    */
    template <class T>
    DiscountedRidgeRegression<T>
    ::DiscountedRidgeRegression(T penalization, T gamma, T p):
        BaseForecaster<T>(), penalization_(penalization),
        gamma_(gamma), p_(p)
    {
        this->parameter_.Reallocate(3);
        this->parameter_(0) = penalization;
        this->parameter_(1) = gamma;
        this->parameter_(2) = p;
    }


    //! Destructor.
    template <class T>
    DiscountedRidgeRegression<T>::~DiscountedRidgeRegression()
    {
    }


    //! Sets the vector of parameters.
    /*!
      \param[in] parameter new vector of parameters with: the penalization,
      the parameter gamma and the parameter p.
    */
    template <class T>
    void DiscountedRidgeRegression<T>
    ::SetParameter(const Vector<T>& parameter)
    {
        if (parameter.GetLength() != 3)
            throw ErrorArgument("DiscountedRidgeRegression::SetParameter"
                                "(const Vector&)",
                                "Discounted ridge regression has three "
                                "parameters, but "
                                + to_str(parameter.GetLength())
                                + " parameters were given.");
        this->parameter_ = parameter;
        penalization_ = parameter(0);
        gamma_ = parameter(1);
        p_ = parameter(2);
    }


    //! Returns the penalization.
    /*!
      \return The penalization.
    */
    template <class T>
    T DiscountedRidgeRegression<T>
    ::GetPenalization() const
    {
        return penalization_;
    }


    //! Sets the penalization.
    /*!
      \param[in] penalization penalization.
    */
    template <class T>
    void DiscountedRidgeRegression<T>
    ::SetPenalization(T penalization)
    {
        penalization_ = penalization;
        this->parameter_(0) = penalization;
    }


    //! Sets the parameter gamma.
    /*!
      \param[in] gamma the parameter gamma.
    */
    template <class T>
    void DiscountedRidgeRegression<T>
    ::SetGamma(T gamma)
    {
        gamma_ = gamma;
        this->parameter_(1) = gamma;
    }


    //! Prepares the aggregation.
    /*! This method should be called before an aggregation method is
      called.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] observation the observational data, indexed by time and
      location.
    */
    template <class T>
    void DiscountedRidgeRegression<T>
    ::Init(const Vector3<T>& ensemble, const Vector2<T>& observation)
    {
        BaseForecaster<T>::Init(ensemble, observation);
    }


    //! Aggregates in the spin-up period at a given time.
    /*! It carries out the sequential aggregation and returns the aggregated
      forecast together with the aggregation weights. This method only carries
      out the aggregation at one time.
      \param[in] t time index at which the aggregation is carried out.
      \param[in] simulation_data the ensemble data to be aggregated indexed by
      member and location.
      \param[in] observation_data the observational data, indexed by location.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] observation the observational data, indexed by time and
      location.
      \param[out] weight aggregation weights, indexed by member and time.
      \param[out] aggregated the aggregated output, indexed by time and
      member.
    */
    template <class T>
    void DiscountedRidgeRegression<T>
    ::AggregateSpinUpPeriod(int t,
                            const Matrix<T>& simulation_data,
                            const Vector<T>& observation_data,
                            const Vector3<T>& ensemble,
                            const Vector2<T>& observation,
                            Matrix<T>& weight,
                            Vector2<T>& aggregated)
    {
        AggregateZero(t, ensemble, weight, aggregated);
    }


    //! Computes aggregation weights at a given time.
    /*! It computes the aggregation weights at one time. This method updates
      the weights at time \a t based on past weights.
      \param[in] t time index at which the aggregation is carried out.
      \param[in] simulation_data the ensemble data to be aggregated indexed by
      member and location.
      \param[in] observation_data the observational data, indexed by location.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] observation the observational data, indexed by time and
      location.
      \param[out] weight past aggregation weights, indexed by member and time.
      \param[out] weight_vector aggregation weights for time \a t.
    */
    template <class T>
    inline void DiscountedRidgeRegression<T>
    ::ComputeWeight(int t,
                    const Matrix<T>& simulation_data,
                    const Vector<T>& observation_data,
                    const Vector3<T>& ensemble,
                    const Vector2<T>& observation,
                    const Matrix<T>& weight,
                    Vector<T>& weight_vector)
    {
        ComputeWeightNonSequential(t, ensemble, observation, weight_vector);
    }


    //! Computes aggregation weights at a given time.
    /*! It computes the aggregation weights at one time. This method computes
      the weights at time \a t without knowledge of the past weights.
      \param[in] t time index at which the aggregation is carried out.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] observation the observational data, indexed by time and
      location.
      \param[out] weight_vector aggregation weights for time \a t.
    */
    template <class T>
    void DiscountedRidgeRegression<T>
    ::ComputeWeightNonSequential(int t,
                                 const Vector3<T>& ensemble,
                                 const Vector2<T>& observation,
                                 Vector<T>& weight_vector)
    {
        int h, l, i, j;

        int Ns = ensemble.GetLength();

        T beta;
        Matrix<T> simulation_data_tmp;
        Vector<T> observation_data_tmp;
        Vector<T> simulation_data_location(Ns);

        Matrix<T, General, RowSymPacked> A(Ns);
        A.SetIdentity();
        Mlt(penalization_, A);
        weight_vector.Zero();
        for (h = 0; h < t; h++)
        {
            beta = T(1) + gamma_ / pow(T(t - h), p_);
            Flatten(ensemble, h, h + 1, simulation_data_tmp);
            Flatten(observation, h, h + 1, observation_data_tmp);
            for (l = 0; l < simulation_data_tmp.GetN(); l++)
            {
                GetCol(simulation_data_tmp, l, simulation_data_location);
                Rank1Update(beta, simulation_data_location, A);

                Add(beta * observation_data_tmp(l), simulation_data_location,
                    weight_vector);
            }
        }

        Vector<int> permutation;
        GetLU(A, permutation);
        SolveLU(A, permutation, weight_vector);
    }


}


#endif
