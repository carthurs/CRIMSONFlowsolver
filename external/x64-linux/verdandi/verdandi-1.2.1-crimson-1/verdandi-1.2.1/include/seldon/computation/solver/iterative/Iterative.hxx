// Copyright (C) 2003-2009 Marc Duruflé
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


#ifndef SELDON_FILE_ITERATIVE_HXX

namespace Seldon
{
  //! Base class for preconditioners
  class Preconditioner_Base
  {
  public :

    Preconditioner_Base();

    // solving M z = r
    template<class Matrix1, class Vector1>
    void Solve(const Matrix1& A, const Vector1 & r, Vector1 & z);

    // solving M^t z = r
    template<class Matrix1, class Vector1>
    void TransSolve(const Matrix1& A, const Vector1& r, Vector1 & z);

  };


  //! Class containing parameters for an iterative resolution
  /*!
    Titer is the precision (float or double), the solved
    linear system can be real or complex
  */
  template<class Titer>
  class Iteration
  {
  protected :
    Titer tolerance; //!< stopping criterion
    Titer facteur_reste; //!< inverse of norm of first residual
    int max_iter; //!< maximum number of iterations
    int nb_iter; //!< number of iterations
    int error_code; //!< error code returned by iterative solver
    bool fail_convergence; //!< true if the iterative solver has converged
    //! print level
    /*!
      0 -> no display
      1 -> displays residual after each 100 iterations
      6 -> displays residual after each iteration
    */
    int print_level;
    bool init_guess_null; //!< true if initial guess is null
    int type_solver; //!< iterative solver used
    int parameter_restart; //!< restart parameter (for Gmres and Gcr)
    int type_preconditioning; //!< preconditioner used

  public :

    Iteration();
    Iteration(int max_iteration, const Titer& tol);
    Iteration(const Iteration<Titer>& outer);

    int GetTypeSolver() const;
    int GetRestart() const;
    Titer GetFactor() const;
    Titer GetTolerance() const;
    int GetNumberIteration() const;

    void SetSolver(int type_resolution, int param_restart, int type_prec);
    void SetRestart(int m);
    void SetTolerance(Titer stopping_criterion);
    void SetMaxNumberIteration(int max_iteration);
    void SetNumberIteration(int nb);

    void ShowMessages();
    void ShowFullHistory();
    void HideMessages();

    template<class Vector1>
    int Init(const Vector1& r);
    bool First() const;

    bool IsInitGuess_Null() const;
    void SetInitGuess(bool type) { init_guess_null = type; }

    template<class Vector1>
    bool Finished(const Vector1& r) const;
    bool Finished(const Titer& r) const;

    void Fail(int i, const string& s);

    Iteration& operator ++ (void);

    int ErrorCode() const;

  };

} // end namespace

#define SELDON_FILE_ITERATIVE_HXX
#endif
