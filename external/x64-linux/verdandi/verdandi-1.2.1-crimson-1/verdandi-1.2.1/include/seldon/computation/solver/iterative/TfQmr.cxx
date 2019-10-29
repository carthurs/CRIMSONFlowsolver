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


#ifndef SELDON_FILE_ITERATIVE_TFQMR_CXX

namespace Seldon
{

  //! Solves a linear system by using Transpose Free Quasi-Minimal Residual
  /*!
    Solves the unsymmetric linear system Ax = b using TFQMR.

    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.

    See: R. W. Freund, A Transpose-Free Quasi-Minimal Residual algorithm for
    non-Hermitian linear system. SIAM J. on Sci. Comp. 14(1993), pp. 470-482

    \param[in] A  Complex General Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] iter Iteration parameters
  */
  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int TfQmr(Matrix1& A, Vector1& x, const Vector1& b,
	    Preconditioner& M, Iteration<Titer> & iter)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;

    Vector1 tmp(b), r0(b), v(b), h(b), w(b),
      y1(b), g(b), y0(b), rtilde(b), d(b);

    typedef typename Vector1::value_type Complexe;
    Complexe sigma, alpha, beta, eta, rho, rho0;
    Titer c, kappa, tau, theta;

    tmp.Zero();
    // x is initial value
    // 1. r0 = M^{-1} (b - A x)
    Copy(b, tmp);
    if (!iter.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), tmp);
    else
      x.Zero();

    M.Solve(A, tmp, r0);

    // we initialize iter
    int success_init = iter.Init(r0);
    if (success_init != 0)
      return iter.ErrorCode();

    // 2. w=y=r
    Copy(r0, w);
    Copy(r0, y1);

    // 3. g=v=M^{-1}Ay
    Copy(y1, v);
    Mlt(A, v, tmp);
    M.Solve(A, tmp, v);
    ++iter;

    Copy(v, g);

    // 4. d=0
    d.Zero();

    // 5. tau = ||r||
    tau = Norm2(r0);

    // 6. theta = eta = 0
    theta = 0.0;
    eta = 0.0;

    // 7. rtilde = r
    Copy(r0, rtilde);

    // 8. rho=dot(rtilde, r)
    rho = DotProd(rtilde, r0);
    rho0 = rho;

    iter.SetNumberIteration(0);
    // Loop until the stopping criteria are reached
    for (;;)
      {
	// 9. 10. 11.
	// sigma=dot(rtilde,v)
	// alpha=rho/sigma
	// y2k=y(2k-1)-alpha*v
	sigma = DotProd(rtilde, v);

	if (sigma == Complexe(0))
	  {
	    iter.Fail(1, "Tfqmr breakdown: sigma=0");
	    break;
	  }
	alpha = rho / sigma;

	//y0 = y1 - alpha * v;
	Copy(y1, y0);
	Add(-alpha, v, y0);

	// 12. h=M^{-1}*A*y
	Copy(y0, h);
	Mlt(A, h, tmp);
	M.Solve(A, tmp, h);
	//split the loop of "for m = 2k-1, 2k"

	//The first one
	// 13. w = w-alpha*M^{-1} A y0
	//w = w - alpha * g;
	Add(-alpha, g, w);
	// 18. d=y0+((theta0^2)*eta0/alpha)*d         //need check breakdown
	if (alpha == Complexe(0))
	  {
	    iter.Fail(2, "Tfqmr breakdown: alpha=0");
	    break;
	  }
	//d = y1 + ( theta * theta * eta / alpha ) * d;
	Mlt(theta * theta * eta / alpha,d);
	Add(Complexe(1), y1, d);

	// 14. theta=||w||_2/tau0       //need check breakdown
	if (tau == Titer(0))
	  {
	    iter.Fail(3, "Tfqmr breakdown: tau=0");
	    break;
	  }
	theta  = Norm2(w) / tau;

	// 15. c = 1/sqrt(1+theta^2)
	c = Titer(1) / sqrt(Titer(1) + theta * theta);

	// 16. tau = tau0*theta*c
	tau = tau * c * theta;

	// 17.  eta = (c^2)*alpha
	eta = c * c * alpha;

	// 19. x = x+eta*d
	Add(eta, d, x);
	// 20. kappa = tau*sqrt(m+1)
	kappa = tau * sqrt( 2.* (iter.GetNumberIteration()+1) );

	// 21. check stopping criterion
	if ( iter.Finished(kappa) )
	  {
	    break;
	  }
	++iter;
	// g = h;
	Copy(h, g);
	//The second one

	// 13. w = w-alpha*M^{-1} A y0
	// w = w - alpha * g;
	Add(-alpha,g,w);
	// 18. d = y0+((theta0^2)*eta0/alpha)*d
	Mlt(theta * theta * eta / alpha,d);
	Add(Complexe(1), y0, d);
	// 14. theta=||w||_2/tau0
	if (tau == Titer(0))
	  {
	    iter.Fail(4, "Tfqmr breakdown: tau=0");
	    break;
	  }
	theta = Norm2(w) / tau;

	// 15. c = 1/sqrt(1+theta^2)
	c = Titer(1) / sqrt(Titer(1) + theta * theta);

	// 16. tau = tau0*theta*c
	tau = tau * c * theta;

	// 17.  eta = (c^2)*alpha
	eta = c * c * alpha;

	// 19. x = x+eta*d
	Add(eta, d, x);

	// 20. kappa = tau*sqrt(m+1)
	kappa = tau * sqrt(2.* (iter.GetNumberIteration()+1)  + 1.);

	// 21. check stopping criterion
	if ( iter.Finished(kappa) )
	  {
	    break;
	  }

	// 22. rho = dot(rtilde,w)
	// 23. beta = rho/rho0                     //need check breakdown

	rho0 = rho;
	rho = DotProd(rtilde, w);
	if (rho0 == Complexe(0))
	  {
	    iter.Fail(5, "tfqmr breakdown: beta=0");
	    break;
	  }
	beta = rho/rho0;

	// 24. y = w+beta*y0
	Copy(w, y1);
	Add(beta, y0, y1);

	// 25. g=M^{-1} A y
	Copy(y1, g);
	Mlt(A, g, tmp);
	M.Solve(A, tmp, g);

	// 26. v = M^{-1}A y + beta*( M^{-1} A y0 + beta*v)

	Mlt(beta*beta, v);
	Add(beta, h, v);
	Add(Complexe(1), g, v);

	++iter;
      }

    return iter.ErrorCode();
  }

} // end namespace

#define SELDON_FILE_ITERATIVE_TFQMR_CXX
#endif
