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


#ifndef SELDON_FILE_ITERATIVE_LSQR_CXX

namespace Seldon
{

  //! Solves a linear system by using Least Squares (LSQR)
  /*!
    Solves the unsymmetric linear system A x = b.

    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.

    \param[in] A  Complex General Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] iter Iteration parameters
  */
  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int Lsqr(Matrix1& A, Vector1& x, const Vector1& b,
	   Preconditioner& M, Iteration<Titer> & iter)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;

    typedef typename Vector1::value_type Complexe;
    Complexe rho, rho_bar, phi, phi_bar, theta, c, s, tmp;
    Titer beta, alpha, rnorm;
    Vector1 v(b), v1(b), u(b), u1(b), w(b);

    int success_init = iter.Init(b);
    if (success_init != 0)
      return iter.ErrorCode();

    Copy(b, u);
    if (!iter.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), u);
    else
      x.Zero();

    rnorm = Norm2(u);

    Copy(b, u);
    beta = Norm2(u);
    tmp = 1.0/beta; Mlt(tmp, u);
    // matrix vector product
    Mlt(SeldonTrans,A, u, v);
    alpha = Norm2(v);
    tmp = 1.0/alpha; Mlt(tmp, v);

    Copy(v,w); x.Zero();

    phi_bar = beta; rho_bar = alpha;

    iter.SetNumberIteration(0);
    // Loop until the stopping criteria are satisfied
    while (! iter.Finished(rnorm))
      {
	// matrix vector product u1 = A*v
	Mlt(A, v, u1);
	// u1 = u1 - alpha*u
	Add(-alpha, u, u1);
	beta = Norm2(u1);
	if (beta == Complexe(0) )
	  {
	    iter.Fail(1, "Lsqr breakdown #1");
	    break;
	  }
	tmp = 1.0/beta; Mlt(tmp, u1);

	// matrix vector  product v1 = A^t u1
	Mlt(SeldonTrans, A, u1, v1);
	// v1 = v1 - beta*v
	Add(-beta, v, v1);
	alpha = Norm2(v1);
	if (alpha == Complexe(0) )
	  {
	    iter.Fail(2, "Lsqr breakdown #2");
	    break;
	  }
	tmp = 1.0/alpha; Mlt(tmp, v1);

	rho = sqrt(rho_bar*rho_bar+beta*beta);
	if (rho == Complexe(0) )
	  {
	    iter.Fail(3, "Lsqr breakdown #3");
	    break;
	  }

	c       = rho_bar/rho;
	s       = beta / rho;
	theta   = s*alpha;
	rho_bar = -c * alpha;
	phi     = c * phi_bar;
	phi_bar = s * phi_bar;

	// x = x + (phi/rho) w
	tmp = phi/rho;
	Add(tmp, w, x);
	// w = v1 - (theta/rho) w
	tmp  = -theta/rho;
	Mlt(tmp,w); Add(Complexe(1), v1, w);

	rnorm = abs(phi_bar);

	Swap(u1,u); Swap(v1,v);

	++iter;
	++iter;
      }

    return iter.ErrorCode();
  }

} // end namespace

#define SELDON_FILE_ITERATIVE_LSQR_CXX
#endif
