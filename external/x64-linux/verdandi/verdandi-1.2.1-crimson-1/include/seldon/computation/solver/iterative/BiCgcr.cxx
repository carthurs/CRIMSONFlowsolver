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


#ifndef SELDON_FILE_ITERATIVE_BICGCR_CXX

namespace Seldon
{

  //! Solves a linear system by using BiCgCr
  /*!
    Solves the symmetric linear system Ax = b using the
    BiConjugate Gradient Conjugate Residual method.

    See Iterative Methods for the Solution of Very Large Complex Symmetric Linear Systems of
    Equations in Electrodynamics,
    Markus Clemens, Thomas Weiland
    Fachbereich 18 Elektrische Nachrichtentechnik, Fachgebiet Theorie Elektromagnetischer Felder,
    Technische Hochschule Darmstadt, Schlossgartenstr. 8, D-64289 Darmstadt, Germany

    \param[in] A  Complex Symmetric Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] iter Iteration parameters
  */
  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int BiCgcr(Matrix1& A, Vector1& x, const Vector1& b,
	     Preconditioner& M, Iteration<Titer> & iter)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;

    typedef typename Vector1::value_type Complexe;
    Complexe rho, mu, alpha, beta, tau;
    Vector1 v(b), w(b), s(b), z(b), p(b), a(b);
    v.Zero(); w.Zero(); s.Zero(); z.Zero();  p.Zero(); a.Zero();

    // we initialize iter
    int success_init = iter.Init(b);
    if (success_init != 0)
      return iter.ErrorCode();

    // we compute the residual v = b - Ax
    Copy(b, v);
    if (!iter.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), v);
    else
      x.Zero();

    iter.SetNumberIteration(0);
    // s = M*v   p = s
    M.Solve(A, v, s); Copy(s, p);
    // a = A*p   w = M*a
    Mlt(A, p, a); M.Solve(A, a, w);
    // we made one product matrix vector
    ++iter;

    while (! iter.Finished(v))
      {
	rho = DotProd(w, v);
	mu = DotProd(w, a);

	// new iterate x = x + alpha*p0  where alpha = (r1,r0)/(p1,p1)
	if (mu == Complexe(0))
	  {
	    iter.Fail(1, "Bicgcr breakdown #1");
	    break;
	  }
	alpha = rho/mu;
	Add(alpha, p, x);

	// new residual r0 = r0 - alpha * p1
	// r1 = r1 - alpha*p2
	Add(-alpha, a, v);
	Add(-alpha, w, s);

	Mlt(A, s, z);
	tau = DotProd(w, z);

	if (tau == Complexe(0))
	  {
	    iter.Fail(2, "Bicgcr breakdown #2");
	    break;
	  }

	beta = tau/mu;

	Mlt(Complexe(-beta), p);
	Add(Complexe(1), s, p);
	Mlt(Complexe(-beta), a);
	Add(Complexe(1), z, a);

	M.Solve(A, a, w);

	++iter;
      }

    return iter.ErrorCode();
  }

} // end namespace

#define SELDON_FILE_ITERATIVE_BICGCR_CXX
#endif
