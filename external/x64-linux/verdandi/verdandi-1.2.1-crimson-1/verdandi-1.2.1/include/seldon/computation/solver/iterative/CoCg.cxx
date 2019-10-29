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


#ifndef SELDON_FILE_ITERATIVE_COCG_CXX

namespace Seldon
{

  //! Solves a linear system by using Conjugate Orthogonal Conjugate Gradient
  /*!
    Solves the symmetric complex linear system A x = b.

    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.

    See H. Van der Vorst, J. Melissen,
    A Petrow-Galerkin type method solving Ax=b where A is symmetric complex
    IEEE Trans. Mag., vol 26, no 2, pp 706-708, 1990

    \param[in] A  Complex Symmetric Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Right hand side of the linear system
    \param[in] M Left preconditioner
    \param[in] iter Iteration parameters
  */
  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int CoCg(Matrix1& A, Vector1& x, const Vector1& b,
	   Preconditioner& M, Iteration<Titer> & iter)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;

    typedef typename Vector1::value_type Complexe;
    Complexe rho, rho_1(0), alpha, beta, delta, zero;
    zero = b(0)*Titer(0);
    rho = zero+Titer(1);

    Vector1 p(b), q(b), r(b), z(b);
    p.Fill(zero); q.Fill(zero); r.Fill(zero); z.Fill(zero);

    // for implementation see Cg
    // we initialize iter
    int success_init = iter.Init(b);
    if (success_init != 0)
      return iter.ErrorCode();

    Copy(b,r);
    if (!iter.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), r);
    else
      x.Fill(zero);

    iter.SetNumberIteration(0);
    // Loop until the stopping criteria are reached
    while (! iter.Finished(r))
      {
	// preconditioning
	M.Solve(A, r, z);

	// instead of (bar(r),z) in CG we compute (r,z)
	rho = DotProd(r, z);

	if (rho == zero)
	  {
	    iter.Fail(1, "Cocg breakdown #1");
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
	// product matrix vector
	Mlt(A, p, q);

	delta = DotProd(p, q);
	if (delta == zero)
	  {
	    iter.Fail(2, "Cocg breakdown #2");
	    break;
	  }
	alpha = rho / delta;

	Add(alpha, p, x);
	Add(-alpha, q, r);

	rho_1 = rho;

	++iter;
      }

    return iter.ErrorCode();
  }


} // end namespace

#define ITERATIVE_COCG_CXX
#endif
