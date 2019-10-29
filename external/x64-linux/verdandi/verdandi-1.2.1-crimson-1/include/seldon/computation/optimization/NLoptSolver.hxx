// Copyright (C) 2011 INRIA
// Author(s): Marc Fragu
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.


#ifndef SELDON_COMPUTATION_OPTIMIZATION_NLOPTSOLVER_HXX
#define SELDON_COMPUTATION_OPTIMIZATION_NLOPTSOLVER_HXX


#include "NLopt.hxx"


namespace Seldon
{


  /////////////////
  // NLOPTSOLVER //
  /////////////////


  //! NLopt optimization.
  class NLoptSolver
  {

  protected:

    typedef double (*cost_ptr)(const Vector<double>&,
                               Vector<double>&, void*);

    //! NLopt optimization solver.
    nlopt::SeldonOpt opt_;
    //! Optimization algorithm.
    nlopt::algorithm algorithm_;
    //! Relative tolerance on the optimization parameters.
    double parameter_tolerance_;
    //! Relative tolerance on the cost function.
    double cost_function_tolerance_;
    /*! \brief Maximum number of function evaluations. It is ignored if it is
      non-positive. */
    int Niteration_max_;
    /*! \brief The vector that stores parameters values. Before optimization,
      stores the initial parameter vector; after optimization, it returns the
      optimized parameters. */
    Vector<double> parameter_;
    /*! \brief The vector that stores gradient values. Before optimization,
      unspecified; after optimization, it returns the gradient vector for
      optimized parameters. */
    Vector<double> gradient_;
    //! The value of cost function for given parameter values.
    double cost_;

  public:
    // Constructor and destructor.
    NLoptSolver();
    ~NLoptSolver();

    void Initialize(int Nparameter, string algorithm,
                    double parameter_tolerance = 1.e-6,
                    double cost_function_tolerance = 1.e-6,
                    int Niteration_max = -1);
    void SetLowerBound(const Vector<double>&);
    void SetUpperBound(const Vector<double>&);
    void SetParameterTolerance(double);
    void SetCostFunctionTolerance(double);
    void SetNiterationMax(int);
    void GetParameterTolerance(double&) const;
    void GetCostFunctionTolerance(double&) const;
    void GetNiterationMax(int&) const;
    void SetParameter(const Vector<double>& parameter);
    void GetParameter(Vector<double>& parameter) const;
    void Optimize(cost_ptr cost, void* argument);
    double GetCost() const;

  };


} // namespace Seldon.


#endif
