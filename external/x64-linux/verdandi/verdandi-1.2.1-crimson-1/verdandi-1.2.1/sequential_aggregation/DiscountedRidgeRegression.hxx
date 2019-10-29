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


#ifndef VERDANDI_FILE_SEQUENTIAL_AGGREGATION_DISCOUNTEDRIDGEREGRESSION_HXX
#define VERDANDI_FILE_SEQUENTIAL_AGGREGATION_DISCOUNTEDRIDGEREGRESSION_HXX


#include "BaseForecaster.hxx"


namespace Verdandi
{


    template <class T>
    class DiscountedRidgeRegression: public BaseForecaster<T>
    {
    protected:
        //! Penalization.
        T penalization_;
        //! Parameter of the sequence of discount factors.
        T gamma_;
        //! Parameter of the sequence of discount factors.
        T p_;

    public:
        DiscountedRidgeRegression();
        DiscountedRidgeRegression(T penalization, T gamma, T p);
        virtual ~DiscountedRidgeRegression();

        virtual void SetParameter(const Vector<T>& parameter);
        T GetPenalization() const;
        void SetPenalization(T penalization);
        void SetGamma(T gamma);

        using BaseForecaster<T>::Aggregate;
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
    };


}


#endif
