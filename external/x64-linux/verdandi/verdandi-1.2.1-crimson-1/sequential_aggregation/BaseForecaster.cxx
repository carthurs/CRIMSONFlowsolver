// Copyright (C) 2008-2011 INRIA
// Author(s): Sébastien Gerchinovitz, Vivien Mallet
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


#ifndef VERDANDI_FILE_SEQUENTIAL_AGGREGATION_BASEFORECASTER_CXX
#define VERDANDI_FILE_SEQUENTIAL_AGGREGATION_BASEFORECASTER_CXX


#include "BaseForecaster.hxx"

#include "seldon/vector/Vector2.cxx"
#include "seldon/vector/Vector3.cxx"
#include "Aggregate.cxx"


namespace Verdandi
{


    //! Main constructor.
    /*! The width of the learning time window is set to 1. */
    template <class T>
    BaseForecaster<T>::BaseForecaster():
        Nspinup_(1), is_convex_(true)
    {
    }


    //! Destructor.
    template <class T>
    BaseForecaster<T>::~BaseForecaster()
    {
    }


    //! Returns the width of the spin-up period.
    /*!
      \return The width of the spin-up period.
    */
    template <class T>
    inline int BaseForecaster<T>::GetNspinup() const
    {
        return Nspinup_;
    }


    //! Sets the width of the spin-up period.
    /*!
      \param[in] Nspinup new width of the spin-up period.
    */
    template <class T>
    void BaseForecaster<T>::SetNspinup(int Nspinup)
    {
        Nspinup_ = Nspinup;
    }


    //! Return the number of parameters.
    /*!
      \return The number of parameters in the algorithm.
    */
    template <class T>
    inline int BaseForecaster<T>::GetNparameter() const
    {
        return parameter_.GetLength();
    }


    //! Returns the vector of parameters.
    /*!
      \return The vector of parameters.
    */
    template <class T>
    inline const Vector<T>&
    BaseForecaster<T>::GetParameter() const
    {
        return parameter_;
    }


    //! Sets the vector of parameters.
    /*!
      \param[in] parameter new vector of parameters.
    */
    template <class T>
    void BaseForecaster<T>::SetParameter(const Vector<T>& parameter)
    {
        parameter_ = parameter;
    }


    //! Is the method carrying out convex aggregation?
    /*! It returns true if the methods produces convex combinations.
      \return true if the method is convex, false otherwise.
    */
    template <class T>
    inline bool BaseForecaster<T>::IsConvex() const
    {
        return is_convex_;
    }


    //! Aggregates the ensemble data.
    /*! It carries out the sequential aggregation and returns the aggregated
      forecast.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] observation the observational data, indexed by time and
      location.
      \param[out] aggregated the aggregated output, indexed by time and
      location.
    */
    template <class T>
    void BaseForecaster<T>
    ::Aggregate(const Vector3<T>& ensemble,
                const Vector2<T>& observation,
                Vector2<T>& aggregated)
    {
        Matrix<T> weight;
        Aggregate(ensemble, observation, weight, aggregated);
    }


    //! Aggregates the ensemble data.
    /*! It carries out the sequential aggregation with a threshold on the
      weights and returns the aggregated forecast.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] observation the observational data, indexed by time and
      location.
      \param[in] threshold_value threshold on the weights.
      \param[out] aggregated the aggregated output, indexed by time and
      location.
    */
    template <class T>
    void BaseForecaster<T>
    ::AggregateThreshold(const Vector3<T>& ensemble,
                         const Vector2<T>& observation,
                         T threshold_value,
                         Vector2<T>& aggregated)
    {
        Matrix<T> weight;
        AggregateThreshold(ensemble, observation,
                           threshold_value, weight, aggregated);
    }


    //! Aggregates the ensemble data.
    /*! It carries out the sequential aggregation and returns the aggregated
      forecast together with the aggregation weights.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] observation the observational data, indexed by time and
      location.
      \param[out] weight aggregation weights, indexed by member and time.
      \param[out] aggregated the aggregated output, indexed by time and
      member.
    */
    template <class T>
    void BaseForecaster<T>
    ::Aggregate(const Vector3<T>& ensemble,
                const Vector2<T>& observation,
                Matrix<T>& weight,
                Vector2<T>& aggregated)
    {
        this->Init(ensemble, observation);
        aggregated = observation;

        int Ns = ensemble.GetLength();
        int Nt = observation.GetLength();

        Matrix<T> simulation_data;
        Vector<T> observation_data;
        weight.Reallocate(Nt, Ns);

        // In the first steps.
        if (Nt != 0)
            this->AggregateSpinUpPeriod(0, simulation_data, observation_data,
                                        ensemble, observation,
                                        weight, aggregated);

        for (int t = 1; t < Nspinup_; t++)
        {
            Flatten(ensemble, 0, t, simulation_data);
            Flatten(observation, 0, t, observation_data);
            this->AggregateSpinUpPeriod(t, simulation_data, observation_data,
                                        ensemble, observation,
                                        weight, aggregated);
        }

        // After these steps.
        for (int t = Nspinup_; t < observation.GetLength(); t++)
        {
            Flatten(ensemble, 0, t, simulation_data);
            Flatten(observation, 0, t, observation_data);
            this->Aggregate(t, simulation_data, observation_data,
                            ensemble, observation, weight, aggregated);
        }
    }


    //! Aggregates the ensemble data.
    /*! It carries out the sequential aggregation with a threshold on the
      weights and returns the aggregated forecast together with the
      aggregation weights.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] observation the observational data, indexed by time and
      location.
      \param[in] threshold_value threshold on the weights.
      \param[out] weight aggregation weights, indexed by member and time.
      \param[out] aggregated the aggregated output, indexed by time and
      member.
    */
    template <class T>
    void BaseForecaster<T>
    ::AggregateThreshold(const Vector3<T>& ensemble,
                         const Vector2<T>& observation,
                         T threshold_value,
                         Matrix<T>& weight,
                         Vector2<T>& aggregated)
    {
        this->Init(ensemble, observation);
        if (is_convex_
            && threshold_value > T(1) / T(ensemble.GetLength()))
            throw ErrorProcessing("AggregatingMethod::AggregateThreshold",
                                  string("Aggregation procedure is convex ")
                                  + "but threshold value is greater than 1 / "
                                  + to_str(ensemble.GetLength()) + ".");

        aggregated = observation;

        int Ns = ensemble.GetLength();
        int Nt = observation.GetLength();

        Matrix<T> simulation_data;
        Vector<T> observation_data;
        weight.Reallocate(Nt, Ns);
        Matrix<T> base_weight(Nt, Ns);

        // In the first steps.
        if (Nt != 0)
            this->AggregateSpinUpPeriodThreshold(0, simulation_data,
                                                 observation_data,
                                                 ensemble, observation,
                                                 threshold_value,
                                                 weight, base_weight,
                                                 aggregated);
        for (int t = 1; t < Nspinup_; t++)
        {
            Flatten(ensemble, 0, t, simulation_data);
            Flatten(observation, 0, t, observation_data);
            this->AggregateSpinUpPeriodThreshold(t, simulation_data,
                                                 observation_data,
                                                 ensemble, observation,
                                                 threshold_value,
                                                 weight, base_weight,
                                                 aggregated);
        }
        // After these steps.
        for (int t = Nspinup_; t < observation.GetLength(); t++)
        {
            Flatten(ensemble, 0, t, simulation_data);
            Flatten(observation, 0, t, observation_data);
            this->AggregateThreshold(t, simulation_data, observation_data,
                                     ensemble, observation,
                                     threshold_value, weight,
                                     base_weight, aggregated);
        }
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
    void BaseForecaster<T>
    ::Init(const Vector3<T>& ensemble, const Vector2<T>& observation)
    {
        int Ns = ensemble.GetLength();
        if (Ns == 0)
            throw ErrorProcessing("BaseForecaster::Init",
                                  "The ensemble has no member to aggregate.");
        for (int s = 0; s < Ns; s++)
            if (ensemble.GetLength(s) != observation.GetLength())
                throw ErrorProcessing("BaseForecaster::Init",
                                      "There is not an equal number "
                                      "of time steps in the " + to_str(s)
                                      + "-th simulation and in the "
                                      "observations.");
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
    void BaseForecaster<T>
    ::AggregateSpinUpPeriod(int t,
                            const Matrix<T>& simulation_data,
                            const Vector<T>& observation_data,
                            const Vector3<T>& ensemble,
                            const Vector2<T>& observation,
                            Matrix<T>& weight,
                            Vector2<T>& aggregated)
    {
        AggregateMean(t, simulation_data, ensemble, weight, aggregated);
    }


    //! Aggregates in the spin-up period at a given time.
    /*! It carries out the sequential aggregation with a threshold on the
      weights and returns the aggregated forecast together with the
      aggregation weights. This method only carries out the aggregation at one
      time.
      \param[in] t time index at which the aggregation is carried out.
      \param[in] simulation_data the ensemble data to be aggregated indexed by
      member and location.
      \param[in] observation_data the observational data, indexed by location.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] observation the observational data, indexed by time and
      location.
      \param[in] threshold_value threshold on the weights.
      \param[out] weight aggregation weights, indexed by member and time.
      \param[out] base_weight aggregation weights before the threshold is
      applied, indexed by member and time.
      \param[out] aggregated the aggregated output, indexed by time and
      member.
    */
    template <class T>
    void BaseForecaster<T>
    ::AggregateSpinUpPeriodThreshold(int t,
                                     const Matrix<T>& simulation_data,
                                     const Vector<T>& observation_data,
                                     const Vector3<T>& ensemble,
                                     const Vector2<T>& observation,
                                     T& threshold_value,
                                     Matrix<T>& weight,
                                     Matrix<T>& base_weight,
                                     Vector2<T>& aggregated)
    {
        this->AggregateSpinUpPeriod(t, simulation_data, observation_data,
                                    ensemble, observation, base_weight,
                                    aggregated);

        int Ns = ensemble.GetLength();

        Vector<T> weight_vector(Ns);
        GetRow(base_weight, t, weight_vector);

        // Thresholding.
        int s;
        if (is_convex_)
        {
            T weight_vector_sum(0);
            for (s = 0; s < Ns; s++)
                if (weight_vector(s) < threshold_value)
                    weight_vector(s) = T(0);
                else
                    weight_vector_sum += weight_vector(s);

            for (s = 0; s < Ns; s++)
                weight_vector(s) = weight_vector(s) / weight_vector_sum;
        }
        else
            for (s = 0; s < Ns; s++)
                if (abs(weight_vector(s)) < threshold_value)
                    weight_vector(s) = T(0);

        SetRow(weight_vector, t, weight);

        // Aggregation.
        ComputeAggregatedValue(t, ensemble, weight_vector, aggregated);
    }


    //! Computes aggregation weights at a given time.
    /*! It computes the aggregation weights at one time. This method updates
      the weights at time \a t based on past weights, which may be more
      efficient than computing the weights from scratch.
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
    void BaseForecaster<T>
    ::ComputeWeight(int t,
                    const Matrix<T>& simulation_data,
                    const Vector<T>& observation_data,
                    const Vector3<T>& ensemble,
                    const Vector2<T>& observation,
                    const Matrix<T>& weight,
                    Vector<T>& weight_vector)
    {
        int Ns = ensemble.GetLength();
        for (int s = 0; s < Ns; s++)
            weight_vector(s) = T(1) / T(Ns);
    }


    //! Computes aggregation weights at a given time.
    /*! It computes the aggregation weights at one time. This method computes
      the weights at time \a t without knowledge of the past weights, which
      may be less efficient than relying on past weights.
      \param[in] t time index at which the aggregation is carried out.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] observation the observational data, indexed by time and
      location.
      \param[out] weight_vector aggregation weights for time \a t.
    */
    template <class T>
    void BaseForecaster<T>
    ::ComputeWeightNonSequential(int t,
                                 const Vector3<T>& ensemble,
                                 const Vector2<T>& observation,
                                 Vector<T>& weight_vector)
    {
        this->Init(ensemble, observation);
        Vector2<T> aggregated = observation;

        int Ns = ensemble.GetLength();
        int Nt = ensemble.GetLength(0);

        Matrix<T> simulation_data;
        Vector<T> observation_data;
        Matrix<T> weight(t + 1, Ns);

        // In the first steps.
        if (Nt != 0)
            this->AggregateSpinUpPeriod(0, simulation_data, observation_data,
                                        ensemble, observation,
                                        weight, aggregated);

        for (int h = 1; h < Nspinup_; h++)
        {
            Flatten(ensemble, 0, h, simulation_data);
            Flatten(observation, 0, h, observation_data);
            this->AggregateSpinUpPeriod(h, simulation_data, observation_data,
                                        ensemble, observation,
                                        weight, aggregated);
        }

        // After these steps.
        for (int h = Nspinup_; h < t + 1; h++)
        {
            Flatten(ensemble, 0, h, simulation_data);
            Flatten(observation, 0, h, observation_data);

            this->ComputeWeight(h, simulation_data, observation_data,
                                ensemble, observation, weight, weight_vector);
            SetRow(weight_vector, h, weight);
        }
    }


    //! Computes the aggregated value based on given weights.
    /*!
      \param[in] t time index at which the aggregation is carried out.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] weight_vector aggregation weights for time \a t.
      \param[out] aggregated the aggregated value, indexed by time and
      location, and updated for time \a t.
    */
    template <class T>
    void BaseForecaster<T>
    ::ComputeAggregatedValue(int t,
                             const Vector3<T>& ensemble,
                             const Vector<T>& weight_vector,
                             Vector2<T>& aggregated)
    {
        int s;
        T aggregated_value;
        for (int l = 0; l < ensemble.GetLength(0, t); l++)
        {
            aggregated_value = T(0);
            for (s = 0; s < ensemble.GetLength(); s++)
                aggregated_value += ensemble(s, t, l) * weight_vector(s);
            aggregated(t, l) = aggregated_value;
        }
    }


    //! Aggregates at a given time.
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
    void BaseForecaster<T>
    ::Aggregate(int t,
                const Matrix<T>& simulation_data,
                const Vector<T>& observation_data,
                const Vector3<T>& ensemble,
                const Vector2<T>& observation,
                Matrix<T>& weight,
                Vector2<T>& aggregated)
    {
        Vector<T> weight_vector(ensemble.GetLength());

        // Weights computation.
        this->ComputeWeight(t, simulation_data, observation_data,
                            ensemble, observation, weight, weight_vector);
        SetRow(weight_vector, t, weight);

        // Aggregation.
        ComputeAggregatedValue(t, ensemble, weight_vector, aggregated);
    }


    //! Aggregates at a given time.
    /*! It carries out the sequential aggregation with a threshold on the
      weights and returns the aggregated forecast together with the
      aggregation weights. This method only carries out the aggregation at one
      time.
      \param[in] t time index at which the aggregation is carried out.
      \param[in] simulation_data the ensemble data to be aggregated indexed by
      member and location.
      \param[in] observation_data the observational data, indexed by location.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] observation the observational data, indexed by time and
      location.
      \param[in] threshold_value threshold on the weights.
      \param[out] weight aggregation weights, indexed by member and time.
      \param[out] base_weight aggregation weights before the threshold is
      applied, indexed by member and time.
      \param[out] aggregated the aggregated output, indexed by time and
      member.
    */
    template <class T>
    void BaseForecaster<T>
    ::AggregateThreshold(int t,
                         const Matrix<T>& simulation_data,
                         const Vector<T>& observation_data,
                         const Vector3<T>& ensemble,
                         const Vector2<T>& observation,
                         T& threshold_value,
                         Matrix<T>& weight,
                         Matrix<T>& base_weight,
                         Vector2<T>& aggregated)
    {
        Vector<T> weight_vector(ensemble.GetLength());

        // Weights computation.
        this->ComputeWeight(t, simulation_data, observation_data,
                            ensemble, observation,
                            base_weight, weight_vector);
        SetRow(weight_vector, t, base_weight);

        // Thresholding.
        int s;
        if (is_convex_)
        {
            T weight_vector_sum(0);
            for (s = 0; s < ensemble.GetLength(); s++)
                if (weight_vector(s) < threshold_value)
                    weight_vector(s) = T(0);
                else
                    weight_vector_sum += weight_vector(s);

            for (s = 0; s < ensemble.GetLength(); s++)
                weight_vector(s) = weight_vector(s) / weight_vector_sum;
        }
        else
            for (s = 0; s < ensemble.GetLength(); s++)
                if (abs(weight_vector(s)) < threshold_value)
                    weight_vector(s) = T(0);

        SetRow(weight_vector, t, weight);

        // Aggregation.
        ComputeAggregatedValue(t, ensemble, weight_vector, aggregated);
    }


    //! Computes the ensemble mean at a given time.
    /*! It carries out the aggregation and returns the aggregated forecast
      together with the aggregation weights. This method only carries out the
      aggregation at one time.
      \param[in] t time index at which the aggregation is carried out.
      \param[in] simulation_data the ensemble data to be aggregated indexed by
      member and location.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[out] weight (equal) aggregation weights, indexed by member and
      time.
      \param[out] aggregated the aggregated output, indexed by time and
      member.
    */
    template <class T>
    void BaseForecaster<T>
    ::AggregateMean(int t,
                    const Matrix<T>& simulation_data,
                    const Vector3<T>& ensemble,
                    Matrix<T>& weight,
                    Vector2<T>& aggregated)
    {
        int s;
        T aggregated_value;
        for (s = 0; s < ensemble.GetLength(); s++)
            weight(t, s) = T(1) / T(ensemble.GetLength());
        for (int l = 0; l < ensemble.GetLength(0, t); l++)
        {
            aggregated_value = T(0);
            for (s = 0; s < ensemble.GetLength(); s++)
                aggregated_value += ensemble(s, t, l);
            aggregated(t, l) = aggregated_value / T(ensemble.GetLength());
        }
    }


    //! Aggregates with zero weights at a given time.
    /*! It carries out the aggregation and returns the aggregated forecast
      together with the aggregation weights. This method only carries out the
      aggregation at one time.
      \param[in] t time index at which the aggregation is carried out.
      \param[in] simulation_data the ensemble data to be aggregated indexed by
      member and location.
      \param[in] observation_data the observational data, indexed by location.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[out] weight aggregation weights (i.e., zeros), indexed by member
      and time.
      \param[out] aggregated the aggregated output, indexed by time and
      member.
    */
    template <class T>
    void BaseForecaster<T>
    ::AggregateZero(int t,
                    const Vector3<T>& ensemble,
                    Matrix<T>& weight,
                    Vector2<T>& aggregated)
    {
        for (int s = 0; s < ensemble.GetLength(); s++)
            weight(t, s) = T(0);
        for (int l = 0; l < ensemble.GetLength(0, t); l++)
            aggregated(t, l) = T(0);
    }


}


#endif
