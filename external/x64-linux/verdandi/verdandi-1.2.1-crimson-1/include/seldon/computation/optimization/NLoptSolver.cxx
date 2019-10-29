// Copyright (C) 2010 INRIA
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


#ifndef SELDON_COMPUTATION_OPTIMIZATION_NLOPTSOLVER_CXX
#define SELDON_COMPUTATION_OPTIMIZATION_NLOPTSOLVER_CXX


#include "NLoptSolver.hxx"


namespace Seldon
{


  //! Default constructor.
  NLoptSolver::NLoptSolver()
  {
  }


  //! Destructor.
  NLoptSolver::~NLoptSolver()
  {
  }


  //! Initializations.
  /*!
    \param[in] Nparameter total number of parameters to be optimized.
    \param[in] algorithm name of the optimization algorithm, one of:
    GN_DIRECT, GN_DIRECT_L, GN_DIRECT_L_RAND, GN_DIRECT_NOSCAL,
    GN_DIRECT_L_NOSCAL, GN_DIRECT_L_RAND_NOSCAL, GN_ORIG_DIRECT,
    GN_ORIG_DIRECT_L, GD_STOGO, GD_STOGO_RAND, LD_LBFGS_NOCEDAL, LD_LBFGS,
    LN_PRAXIS, LD_VAR1, LD_VAR2, LD_TNEWTON, LD_TNEWTON_RESTART,
    LD_TNEWTON_PRECOND, LD_TNEWTON_PRECOND_RESTART, GN_CRS2_LM, GN_MLSL,
    GD_MLSL, GN_MLSL_LDS, GD_MLSL_LDS, LD_MMA, LN_COBYLA, LN_NEWUOA,
    LN_NEWUOA_BOUND, LN_NELDERMEAD, LN_SBPLX, LN_AUGLAG, LD_AUGLAG,
    LN_AUGLAG_EQ, LD_AUGLAG_EQ, LN_BOBYQA, GN_ISRES, AUGLAG, AUGLAG_EQ,
    G_MLSL, G_MLSL_LDS, LD_SLSQP, NUM_ALGORITHMS.
    \param[in] parameter_tolerance relative tolerance on the parameters.
    When the variation of the parameters, after one step of the algorithm, has
    changed by less than \a parameter_tolerance multiplied by the value of the
    parameters, the optimization is stopped. If you do not want to use a
    particular tolerance termination, you can just set that tolerance to zero
    and it will be ignored.
    \param[in] cost_function_tolerance relative tolerance on the cost
    function. When the variation of the cost function, after one step of the
    algorithm, has changed by less than \a cost_function_tolerance multiplied
    by the value of the cost function, the optimization is stopped. If you
    do not want to use a particular tolerance termination, you can just set
    that tolerance to zero and it will be ignored.
    \param[in] Niteration_max maximum number of cost function evaluations.
    It is ignored if it is non-positive.
  */
  void NLoptSolver::Initialize(int Nparameter, string algorithm,
                               double parameter_tolerance,
                               double cost_function_tolerance,
                               int Niteration_max)
  {
    map<string, nlopt::algorithm> algorithm_map;
    map<string, nlopt::algorithm>::iterator it;
    algorithm_map["GN_DIRECT"] = nlopt::GN_DIRECT;
    algorithm_map["GN_DIRECT_L"] = nlopt::GN_DIRECT_L;
    algorithm_map["GN_DIRECT_L_RAND"] = nlopt::GN_DIRECT_L_RAND;
    algorithm_map["GN_DIRECT_NOSCAL"] = nlopt::GN_DIRECT_NOSCAL;
    algorithm_map["GN_DIRECT_L_NOSCAL"] = nlopt::GN_DIRECT_L_NOSCAL;
    algorithm_map["GN_DIRECT_L_RAND_NOSCAL"]
      = nlopt::GN_DIRECT_L_RAND_NOSCAL;
    algorithm_map["GN_ORIG_DIRECT"] = nlopt::GN_ORIG_DIRECT;
    algorithm_map["GN_ORIG_DIRECT_L"] = nlopt::GN_ORIG_DIRECT_L;
    algorithm_map["GD_STOGO"] = nlopt::GD_STOGO;
    algorithm_map["GD_STOGO_RAND"] = nlopt::GD_STOGO_RAND;
    algorithm_map["LD_LBFGS_NOCEDAL"] = nlopt::LD_LBFGS_NOCEDAL;
    algorithm_map["LD_LBFGS"] = nlopt::LD_LBFGS;
    algorithm_map["LN_PRAXIS"] = nlopt::LN_PRAXIS;
    algorithm_map["LD_VAR1"] = nlopt::LD_VAR1;
    algorithm_map["LD_VAR2"] = nlopt::LD_VAR2;
    algorithm_map["LD_TNEWTON"] = nlopt::LD_TNEWTON;
    algorithm_map["LD_TNEWTON_RESTART"] = nlopt::LD_TNEWTON_RESTART;
    algorithm_map["LD_TNEWTON_PRECOND"] = nlopt::LD_TNEWTON_PRECOND;
    algorithm_map["LD_TNEWTON_PRECOND_RESTART"]
      = nlopt::LD_TNEWTON_PRECOND_RESTART;
    algorithm_map["GN_CRS2_LM"] = nlopt::GN_CRS2_LM;
    algorithm_map["GN_MLSL"] = nlopt::GN_MLSL;
    algorithm_map["GD_MLSL"] = nlopt::GD_MLSL;
    algorithm_map["GN_MLSL_LDS"] = nlopt::GN_MLSL_LDS;
    algorithm_map["GD_MLSL_LDS"] = nlopt::GD_MLSL_LDS;
    algorithm_map["LD_MMA"] = nlopt::LD_MMA;
    algorithm_map["LN_COBYLA"] = nlopt::LN_COBYLA;
    algorithm_map["LN_NEWUOA"] = nlopt::LN_NEWUOA;
    algorithm_map["LN_NEWUOA_BOUND"] = nlopt::LN_NEWUOA_BOUND;
    algorithm_map["LN_NELDERMEAD"] = nlopt::LN_NELDERMEAD;
    algorithm_map["LN_SBPLX"] = nlopt::LN_SBPLX;
    algorithm_map["LN_AUGLAG"] = nlopt::LN_AUGLAG;
    algorithm_map["LD_AUGLAG"] = nlopt::LD_AUGLAG;
    algorithm_map["LN_AUGLAG_EQ"] = nlopt::LN_AUGLAG_EQ;
    algorithm_map["LD_AUGLAG_EQ"] = nlopt::LD_AUGLAG_EQ;
    algorithm_map["LN_BOBYQA"] = nlopt::LN_BOBYQA;
    algorithm_map["GN_ISRES"] = nlopt::GN_ISRES;
    algorithm_map["AUGLAG"] = nlopt::AUGLAG;
    algorithm_map["AUGLAG_EQ"] = nlopt::AUGLAG_EQ;
    algorithm_map["G_MLSL"] = nlopt::G_MLSL;
    algorithm_map["G_MLSL_LDS"] = nlopt::G_MLSL_LDS;
    algorithm_map["LD_SLSQP"] = nlopt::LD_SLSQP;

    it = algorithm_map.find(algorithm);
    if (it == algorithm_map.end())
      WrongArgument("void NLoptSolver::Initialize(int, string, double)",
                    "Unknown algorithm. Implemented algorithms are:"
                    " GN_DIRECT, "
                    "GN_DIRECT_L, GN_DIRECT_L_RAND, "
                    "GN_DIRECT_NOSCAL, "
                    "GN_DIRECT_L_NOSCAL, GN_DIRECT_L_RAND_NOSCAL,"
                    " GN_ORIG_DIRECT, "
                    "GN_ORIG_DIRECT_L, GD_STOGO, GD_STOGO_RAND, "
                    "LD_LBFGS_NOCEDAL, LD_LBFGS, LN_PRAXIS, "
                    "LD_VAR1, LD_VAR2, LD_TNEWTON, "
                    "LD_TNEWTON_RESTART, LD_TNEWTON_PRECOND, "
                    "LD_TNEWTON_PRECOND_RESTART, GN_CRS2_LM, "
                    "GN_MLSL,"
                    "GD_MLSL, GN_MLSL_LDS, GD_MLSL_LDS, "
                    "LD_MMA, LN_COBYLA, LN_NEWUOA, "
                    "LN_NEWUOA_BOUND, LN_NELDERMEAD, "
                    "LN_SBPLX, LN_AUGLAG, LD_AUGLAG, "
                    "LN_AUGLAG_EQ, LD_AUGLAG_EQ, LN_BOBYQA, "
                    "GN_ISRES, AUGLAG, AUGLAG_EQ,"
                    "G_MLSL, G_MLSL_LDS, LD_SLSQP, NUM_ALGORITHMS");
    else
      algorithm_ = it->second;

    opt_ = nlopt::SeldonOpt(algorithm_, Nparameter);

    parameter_tolerance_ = parameter_tolerance;

    cost_function_tolerance_ = cost_function_tolerance;

    Niteration_max_ = Niteration_max;
  }


  //! Sets lower bounds on the parameters.
  /*!
    \param[in] lower_bound the lower bound vector.
  */
  void NLoptSolver::SetLowerBound(const Vector<double>& lower_bound)
  {
    if (lower_bound.GetSize() != 0)
      opt_.set_lower_bounds(lower_bound);
  }


  //! Sets upper bounds on the parameters.
  /*!
    \param[in] upper_bound the lower bound vector.
  */
  void NLoptSolver::SetUpperBound(const Vector<double>& upper_bound)
  {
    if (upper_bound.GetSize() != 0)
      opt_.set_upper_bounds(upper_bound);
  }


  //! Sets the relative tolerance on the parameters.
  /*!
    \param[in] tolerance relative tolerance on the parameters. When the
    variation of every parameter, after one step of the algorithm, has
    changed by less than \a tolerance multiplied by the value of the
    parameter, the optimization is stopped. If you do not want to use a
    particular tolerance termination, you can just set that tolerance to zero
    and it will be ignored.
  */
  void NLoptSolver::SetParameterTolerance(double tolerance)
  {
    parameter_tolerance_ = tolerance;
  }


  //! Sets the relative tolerance on the cost function.
  /*!
    \param[in] tolerance relative tolerance on the cost function. When the
    variation of the cost function, after one step of the algorithm, has
    changed by less than \a tolerance multiplied by the value of the cost
    function, the optimization is stopped. If you do not want to use a
    particular tolerance termination, you can just set that tolerance to zero
    and it will be ignored.
  */
  void NLoptSolver::SetCostFunctionTolerance(double tolerance)
  {
    cost_function_tolerance_ = tolerance;
  }


  //! Sets the maximum number of cost function evaluations.
  /*!
    \param[in] Niteration_max maximum number of cost function evaluations.
    It is ignored if it is non-positive.
  */
  void NLoptSolver::SetNiterationMax(int Niteration_max)
  {
    Niteration_max_ = Niteration_max;
  }


  //! Gets the relative tolerance on the parameters.
  /*!
    \param[out] tolerance relative tolerance on the parameters . When the
    variation of every parameter, after one step of the algorithm, has
    changed by less than \a tolerance multiplied by the value of the
    parameter, the optimization is stopped.
  */
  void NLoptSolver::GetParameterTolerance(double& tolerance) const
  {
    tolerance = parameter_tolerance_;
  }


  //! Gets the relative tolerance on the cost function.
  /*!
    \param[out] tolerance relative tolerance on the cost function . When the
    variation of the cost function, after one step of the algorithm, has
    changed by less than \a tolerance multiplied by the value of the cost
    function, the optimization is stopped.
  */
  void NLoptSolver::GetCostFunctionTolerance(double& tolerance) const
  {
    tolerance = cost_function_tolerance_;
  }


  //! Gets the maximum number of cost function evaluations.
  /*!
    \param[out] Niteration_max maximum number of cost function evaluations.
  */
  void NLoptSolver::GetNiterationMax(int& Niteration_max) const
  {
    Niteration_max = Niteration_max_;
  }


  //! Sets the parameters.
  /*!
    \param[in] parameter the parameters vector.
  */
  void NLoptSolver::SetParameter(const Vector<double>& parameter)
  {
    parameter_.Reallocate(parameter.GetM());
    Copy(parameter, parameter_);
  }


  //! Gets the parameters.
  /*!
    \param[out] parameter the parameters vector.
  */
  void NLoptSolver::GetParameter(Vector<double>& parameter) const
  {
    parameter.Reallocate(parameter_.GetM());
    Copy(parameter_, parameter);
  }


  //! Optimization.
  /*!
    \param[in] cost pointer to the cost function. This function takes as first
    argument a 'const Vector<double>&' of parameters. The second argument of
    the function is a 'Vector<double>&' which must be, on exit, the gradient
    of the cost function. NLopt will allocate it before the call. Note that,
    in case a derivative-free algorithm is used, this gradient vector is
    empty, and the cost function is not supposed to compute it (it is thus
    recommended to test the length of the vector in the cost function). The
    third argument of the function is \a argument, provided as 'void *'. The
    cost function returns the cost value in 'double'.
    \param[in] argument third argument of the cost function. This argument is
    passed to the cost function after the parameters vector and the gradient.
  */
  void NLoptSolver::Optimize(cost_ptr cost, void* argument)
  {
    int Nparameter;
    if (0 == (Nparameter = parameter_.GetM()))
      throw WrongArgument("NLoptSolver::Optimize()",
                          "The vector of parameters to be optimized"
                          " is empty.");

    opt_.set_min_objective(cost, argument);
    opt_.set_xtol_rel(parameter_tolerance_);
    opt_.set_ftol_rel(cost_function_tolerance_);
    opt_.set_maxeval(Niteration_max_);
    nlopt::result result = opt_.optimize(parameter_, cost_);

    if (result < 0)
      throw Error("NLoptSolver::Optimize()", "Nlopt failed.");
  }


  //! Returns the value of the cost function.
  /*! This method should be called after the optimization, and it returns the
    value of the cost function associated with the optimized parameters.
    \return the value of the cost function.
  */
  double NLoptSolver::GetCost() const
  {
    return cost_;
  }


} // namespace Seldon.


#endif
