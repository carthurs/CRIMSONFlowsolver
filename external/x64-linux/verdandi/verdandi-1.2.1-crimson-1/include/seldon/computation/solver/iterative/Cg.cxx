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


#ifndef SELDON_FILE_ITERATIVE_CG_CXX

namespace Seldon
{

  //! Solves a linear system by using Conjugate Gradient (CG)
  /*!
    Solves the symmetric positive definite linear system A x = b.

    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.

    See M. R. Hestenes nd E. Stiefel, Methods of conjugate gradients for
    solving linear system, Journal of Research of the National Bureau of
    Standards, 49(1952), pp. 409-436

    \param[in] A  Real Symmetric Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] iter Iteration parameters
  */
  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int Cg(Matrix1& A, Vector1& x, const Vector1& b,
	 Preconditioner& M, Iteration<Titer> & iter)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;

    typedef typename Vector1::value_type Complexe;
    Complexe rho(1), rho_1(1), alpha, beta,delta;
    Vector1 p(b), q(b), r(b), z(b);

    // we initialize iter
    int success_init = iter.Init(b);
    if (success_init != 0)
      return iter.ErrorCode();

    // we compute the initial residual r = b - Ax
    Copy(b, r);
    if (!iter.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), r);
    else
      x.Zero();

    iter.SetNumberIteration(0);
    // Loop until the stopping criteria are satisfied
    while (! iter.Finished(r))
      {

	// Preconditioning z = M^{-1} r
	M.Solve(A, r, z);

	// rho = (conj(r),z)
	rho = DotProdConj(r, z);

	if (rho == Complexe(0) )
	  {
	    iter.Fail(1, "Cg breakdown #1");
	    break;
	  }

	if (iter.First())
	  Copy(z, p);
	else
	  {
	    // p = beta*p + z  where  beta = rho_i/rho_{i-1}
	    beta = rho / rho_1;
	    Mlt(beta, p);
	    Add(Complexe(1), z, p);
	  }

	// matrix vector product q = A*p
	Mlt(A, p, q);
	delta = DotProdConj(p, q);
	if (delta == Complexe(0))
	  {
	    iter.Fail(2, "Cg breakdown #2");
	    break;
	  }
	alpha = rho / delta;

	// x = x + alpha*p  and r = r - alpha*q  where alpha = rho/(bar(p),q)
	Add(alpha, p, x);
	Add(-alpha, q, r);

	rho_1 = rho;

	++iter;
      }

    return iter.ErrorCode();
  }


} // end namespace

#define SELDON_FILE_ITERATIVE_CG_CXX
#endif

