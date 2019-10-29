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


#ifndef VERDANDI_FILE_SEQUENTIAL_AGGREGATION_BASEFORECASTER_HXX
#define VERDANDI_FILE_SEQUENTIAL_AGGREGATION_BASEFORECASTER_HXX


#include "seldon/vector/Vector2.hxx"
#include "seldon/vector/Vector3.hxx"
#include "Aggregate.hxx"


namespace Verdandi
{


    template <class T>
    class BaseForecaster
    {
    protected:
        //! Number of step in the spin-up period.
        int Nspinup_;

        //! Parameters.
        Vector<T> parameter_;

        //! Is the aggregation procedure convex?
        bool is_convex_;

    public:
        BaseForecaster();
        virtual ~BaseForecaster();

        int GetNspinup() const;
        void SetNspinup(int Nspinup);

        int GetNparameter() const;
        const Vector<T>& GetParameter() const;
        virtual void SetParameter(const Vector<T>& parameter);

        bool IsConvex() const;

        virtual void Aggregate(const Vector3<T>& ensemble,
                               const Vector2<T>& observation,
                               Vector2<T>& aggregated_simulation);
        virtual void
        AggregateThreshold(const Vector3<T>& ensemble,
                           const Vector2<T>& observation,
                           T threshold_value,
                           Vector2<T>& aggregated_simulation);
        virtual void Aggregate(const Vector3<T>& ensemble,
                               const Vector2<T>& observation,
                               Matrix<T>& weight,
                               Vector2<T>& aggregated_simulation);
        virtual void
        AggregateThreshold(const Vector3<T>& ensemble,
                           const Vector2<T>& observation,
                           T threshold_value,
                           Matrix<T>& weight,
                           Vector2<T>& aggregated_simulation);

        virtual void Init(const Vector3<T>& ensemble,
                          const Vector2<T>& observation);

        virtual void
        AggregateSpinUpPeriod(int t,
                              const Matrix<T>& simulation_data,
                              const Vector<T>& observation_data,
                              const Vector3<T>& ensemble,
                              const Vector2<T>& observation,
                              Matrix<T>& weight,
                              Vector2<T>& aggregated_simulation);
        virtual void
        AggregateSpinUpPeriodThreshold(int t,
                                       const Matrix<T>& simulation_data,
                                       const Vector<T>& observation_data,
                                       const Vector3<T>& ensemble,
                                       const Vector2<T>& observation,
                                       T& threshold_value,
                                       Matrix<T>& weight,
                                       Matrix<T>& base_weight,
                                       Vector2<T>& aggregated_simulation);
        virtual void
        ComputeWeight(int t,
                      const Matrix<T>& simulation_data,
                      const Vector<T>& observation_data,
                      const Vector3<T>& ensemble,
                      const Vector2<T>& observation,
                      const Matrix<T>& weight,
                      Vector<T>& weight_vector);
        virtual void
        ComputeWeightNonSequential(int t,
                                   const Vector3<T>& ensemble,
                                   const Vector2<T>& observation,
                                   Vector<T>& weight_vector);
        virtual void
        ComputeAggregatedValue(int t,
                               const Vector3<T>& ensemble,
                               const Vector<T>& weight_vector,
                               Vector2<T>& aggregated_simulation);
        virtual void
        Aggregate(int t,
                  const Matrix<T>& simulation_data,
                  const Vector<T>& observation_data,
                  const Vector3<T>& ensemble,
                  const Vector2<T>& observation,
                  Matrix<T>& weight,
                  Vector2<T>& aggregated_simulation);
        virtual void
        AggregateThreshold(int t,
                           const Matrix<T>& simulation_data,
                           const Vector<T>& observation_data,
                           const Vector3<T>& ensemble,
                           const Vector2<T>& observation,
                           T& threshold_value,
                           Matrix<T>& weight,
                           Matrix<T>& base_weight,
                           Vector2<T>& aggregated_simulation);
        void
        AggregateMean(int t,
                      const Matrix<T>& simulation_data,
                      const Vector3<T>& ensemble,
                      Matrix<T>& weight,
                      Vector2<T>& aggregated_simulation);
        void
        AggregateZero(int t,
                      const Vector3<T>& ensemble,
                      Matrix<T>& weight,
                      Vector2<T>& aggregated_simulation);
    };


}


#endif
