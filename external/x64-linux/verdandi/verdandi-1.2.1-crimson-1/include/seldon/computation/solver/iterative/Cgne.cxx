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


#ifndef SELDON_FILE_ITERATIVE_CGNE_CXX

namespace Seldon
{

  //! Solves a linear system using Conjugate Gradient Normal Equation (CGNE)
  /*!
    Solves the unsymmetric linear system A x = b.

    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.

    \param[in] A  Complex General Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Right hand side of the linear system
    \param[in] M Left preconditioner
    \param[in] iter Iteration parameters
  */
  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int Cgne(Matrix1& A, Vector1& x, const Vector1& b,
	   Preconditioner& M, Iteration<Titer> & iter)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;

    typedef typename Vector1::value_type Complexe;
    Complexe rho(1), rho_1(0), alpha, beta, delta;
    Vector1 p(b), q(b), r(b), z(b);
    Titer dp;

    // x should be equal to 0
    // see Cg to understand implementation
    // we solve A^t A x = A^t b
    // left-preconditioner is equal to M M^t

    // q = A^t b
    Mlt(SeldonTrans, A, b, q);
    // we initialize iter
    int success_init = iter.Init(q);
    if (success_init != 0)
      return iter.ErrorCode();

    if (!iter.IsInitGuess_Null())
      {
	// r = A^t b - A^t A x
	Mlt(A, x, p);
	Mlt(SeldonTrans, A, p, r); Mlt(Complexe(-1), r);
	Add(Complexe(1), q, r);
      }
    else
      {
	Copy(q, r);
	x.Zero();
      }

    dp = Norm2(q);
    iter.SetNumberIteration(0);
    // Loop until the stopping criteria are satisfied
    while (! iter.Finished(dp))
      {
	// Preconditioning
	M.TransSolve(A, r, q);
	M.Solve(A, q, z);

	rho = DotProd(r,z);
	if (rho == Complexe(0))
	  {
	    iter.Fail(1, "Cgne breakdown #1");
	    break;
	  }

	if (iter.First())
	  Copy(z, p);
	else
	  {
	    beta = rho / rho_1;
	    Mlt(beta, p);
	    Add(Complexe(1), z, p);
	  }

	// instead of q = A*p, we compute q = A^t A *p
	Mlt(A, p, q);
	Mlt(SeldonTrans, A, q, z);

	delta = DotProd(p, z);
	if (delta == Complexe(0))
	  {
	    iter.Fail(2, "Cgne breakdown #2");
	    break;
	  }
	alpha = rho / delta;

	Add(alpha, p, x);
	Add(-alpha, z, r);

	rho_1 = rho;
	dp = Norm2(r);

	// two iterations, because of two multiplications with A
	++iter;
	++iter;
      }

    return iter.ErrorCode();
  }


} // end namespace

#define SELDON_FILE_ITERATIVE_CGNE_CXX
#endif
