// Copyright (C) 2003-2009 Marc Duruflé
// Copyright (C) 2001-2009 Vivien Mallet
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


#ifndef SELDON_FILE_ITERATIVE_BICG_CXX

namespace Seldon
{

  //! Solves a linear system by using BiConjugate Gradient (BICG)
  /*!
    Solves the unsymmetric linear system Ax = b
    using the Preconditioned BiConjugate Gradient method.

    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.

    See: R. Fletcher, Conjugate gradient methods for indefinite systems,
    In Numerical Analysis Dundee 1975, G. Watson, ed., Springer Verlag,
    Berlin, New York, 1976 pp. 73-89

    \param[in] A  Complex General Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] iter Iteration parameters
  */
  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int BiCg(Matrix1& A, Vector1& x, const Vector1& b,
	   Preconditioner& M, Iteration<Titer> & iter)
  {
    int N = A.GetM();
    if (N <= 0)
      return 0;

    typedef typename Vector1::value_type Complexe;
    Complexe rho_1, rho_2(0), alpha, beta, delta;

    Vector1 r(b), z(b), p(b), q(b);
    Vector1 r_tilde(b), z_tilde(b), p_tilde(b), q_tilde(b);

    // we initialize iter
    int success_init = iter.Init(b);
    if (success_init != 0)
      return iter.ErrorCode();

    // we compute the residual r = b - Ax
    Copy(b, r);
    if (!iter.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), r);
    else
      x.Zero();

    Copy(r, r_tilde);

    iter.SetNumberIteration(0);
    // Loop until the stopping criteria are reached
    while (! iter.Finished(r))
      {
	// preconditioning z = M^{-1} r and z_tilde = M^{-t} r_tilde
	M.Solve(A, r, z);
	M.TransSolve(A, r_tilde, z_tilde);
	// rho_1 = (z,r_tilde)
	rho_1 = DotProd(z, r_tilde);

	if (rho_1 == Complexe(0))
	  {
	    iter.Fail(1, "Bicg breakdown #1");
	    break;
	  }

	if (iter.First())
	  {
	    Copy(z, p);
	    Copy(z_tilde, p_tilde);
	  }
	else
	  {
	    // p=beta*p+z  where beta = rho_i/rho_{i-1}
	    // p_tilde=beta*p_tilde+z_tilde
	    beta = rho_1 / rho_2;
	    Mlt(beta, p);
	    Add(Complexe(1), z, p);
	    Mlt(beta, p_tilde);
	    Add(Complexe(1), z_tilde, p_tilde);
	  }

	// we do the product matrix vector and transpose matrix vector
	// q = A*p    q_tilde = A^t p_tilde
	Mlt(A, p, q);
	++iter;
	Mlt(SeldonTrans, A, p_tilde, q_tilde);

	delta = DotProd(p_tilde, q);
	if (delta == Complexe(0))
	  {
	    iter.Fail(2, "Bicg breakdown #2");
	    break;
	  }

	alpha = rho_1 / delta;

	// the new iterate x=x+alpha*p and residual r=r-alpha*q
	// where alpha = rho_i/delta
	Add(alpha, p, x);
	Add(-alpha, q, r);
	Add(-alpha, q_tilde, r_tilde);

	rho_2 = rho_1;

	++iter;
      }

    return iter.ErrorCode();
  }

}

#define SELDON_FILE_ITERATIVE_BICG_CXX
#endif
