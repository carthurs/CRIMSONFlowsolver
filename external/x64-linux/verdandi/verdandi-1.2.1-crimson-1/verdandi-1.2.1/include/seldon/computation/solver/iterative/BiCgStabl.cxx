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


#ifndef SELDON_FILE_ITERATIVE_BICGSTABL_CXX

namespace Seldon
{

  //! Implements BiConjugate Gradient Stabilized (BICG-STAB(l))
  /*!
    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.

    See G L.G. Sleijpen, D. R. Fokkema
    BICGSTAB(l) For Linear Equations Involving Unsymmetric Matrices With
    Complex Spectrum

    \param[in] A  Complex General Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] iter Iteration parameters
  */
  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int BiCgStabl(Matrix1& A, Vector1& x, const Vector1& b,
		Preconditioner& M, Iteration<Titer> & iter)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;

    typedef typename Vector1::value_type Complexe;
    int m = iter.GetRestart();
    int l = m;
    Complexe rho_0, rho_1, alpha, beta, omega, sigma, zero, unity;
    zero = 0.0; unity = 1.0;

    // q temporary vector before preconditioning, r0 initial residual
    Vector1 q(b), r0(b);
    Vector<Complexe> gamma(l+1), gamma_prime(l+1), gamma_twice(l+1);
    Matrix<Complexe, General, RowMajor> tau(l+1,l+1);
    // history of u and residual r
    std::vector<Vector1> r(l+1, b), u(l+1, b);
    for (int i = 0; i <= l; i++)
      {
	r[i].Zero();
	u[i].Zero();
      }
    tau.Zero(); gamma.Zero(); gamma_prime.Zero(); gamma_twice.Zero();

    // we compute the residual r = (b - Ax)
    Copy(b, r[0]);
    if (!iter.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), r[0]);
    else
      x.Zero();

    // we initialize iter
    int success_init = iter.Init(b);
    if (success_init !=0 )
      return iter.ErrorCode();

    Copy(r[0], r0); // we keep the first residual

    // we initialize constants
    rho_0 = unity; alpha = zero; omega = unity; tau.Zero();

    iter.SetNumberIteration(0);
    // Loop until the stopping criteria are satisfied
    while (! iter.Finished(r[0]))
      {
	rho_0 *= -omega;

	// Bi-CG Part
	for (int j = 0; j < l; j++)
	  {
	    rho_1 = DotProd(r[j], r0);
	    if (rho_0 == zero)
	      {
		iter.Fail(1, "Bicgstabl breakdown #1");
		break;
	      }
	    beta = alpha*(rho_1/rho_0);
	    rho_0 = rho_1;
	    for (int i = 0; i <= j; i++)
	      {
		Mlt(-beta, u[i]); Add(unity, r[i], u[i]);
	      }
	    M.Solve(A, u[j], q); // preconditioning
	    Mlt(A, q, u[j+1]); // product Matrix Vector

	    ++iter;
	    sigma = DotProd(u[j+1], r0);
	    if (sigma == zero)
	      {
		iter.Fail(2, "Bicgstabl Breakdown #2");
		break;
	      }
	    alpha = rho_1/sigma;
	    Add(alpha, u[0], x);
	    for (int i = 0; i <= j; i++)
	      Add(-alpha, u[i+1], r[i]);

	    M.Solve(A, r[j], q); // preconditioning
	    Mlt(A, q, r[j+1]); // product matrix vector

	    ++iter;
	  }

	// MR Part  modified Gram-Schmidt
	for (int j = 1; j <= l; j++)
	  {
	    for (int i = 1; i < j; i++)
	      {
		if (gamma(i) != zero)
		  {
		    tau(i,j) = DotProd(r[j], r[i])/gamma(i);
		    Add(-tau(i,j), r[i], r[j]);
		  }
	      }
	    gamma(j) = DotProd(r[j], r[j]);
	    if (gamma(j) != zero)
	      gamma_prime(j) = DotProd(r[0], r[j])/gamma(j);
	  }

	// gamma = tau-1 * gamma_prime
	gamma(l) = gamma_prime(l); omega = gamma(l);
	for (int j = l-1; j >= 1; j--)
	  {
	    sigma = zero;
	    for (int i = j+1; i <= l; i++)
	      sigma += tau(j,i)*gamma(i);

	    gamma(j) = gamma_prime(j)-sigma;
	  }

	// gamma_twice=T*S*gamma
	for (int j = 1; j <= l-1; j++)
	  {
	    sigma = zero;
	    for (int i = j+1; i <= l-1; i++)
	      sigma += tau(j,i)*gamma(i+1);

	    gamma_twice(j) = gamma(j+1)+sigma;
	  }

	// update
	Add(gamma(1), r[0], x);
	Add(-gamma_prime(l), r[l], r[0]);
	Add(-gamma(l), u[l], u[0]);
	for (int j = 1;j <= l-1; j++)
	  {
	    Add(-gamma(j), u[j], u[0]);
	    Add(gamma_twice(j), r[j], x);
	    Add(-gamma_prime(j), r[j], r[0]);
	  }
      }
    // change of coordinates (right preconditioning)
    Copy(x,q); M.Solve(A, q, x);
    return iter.ErrorCode();
  }

}

#define SELDON_FILE_ITERATIVE_BICGSTABL_CXX
#endif
