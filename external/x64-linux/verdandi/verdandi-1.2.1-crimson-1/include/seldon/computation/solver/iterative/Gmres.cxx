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


#ifndef SELDON_FILE_ITERATIVE_GMRES_CXX

namespace Seldon
{

  //! Solves a linear system by using Generalized Minimum Residual (GMRES)
  /*!
    Solves the unsymmetric linear system Ax = b using restarted GMRES.

    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.

    See: Y. Saad and M. Schulter. GMRES: A generalized minimum residual
    algorithm for solving nonsysmmetric linear systems, SIAM
    J. Sci. Statist. Comp.  7(1986), pp, 856-869

    \param[in] A  Complex General Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] outer Iteration parameters
  */
  template <class Titer, class MatrixSparse, class Vector1, class Preconditioner>
  int Gmres(MatrixSparse& A, Vector1& x, const Vector1& b,
	    Preconditioner& M, Iteration<Titer> & outer)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;

    typedef typename Vector1::value_type Complexe;
    Complexe zero(0);

    int m = outer.GetRestart();
    // V is the array of orthogonal basis contructed
    // from the Krylov subspace (v0,A*v0,A^2*v0,...,A^m*v0)
    std::vector<Vector1> V(m+1, b);

    // Upper triangular hessenberg matrix
    // we don't store the sub-diagonal
    // we apply rotations to eliminate this sub-diagonal
    Matrix<Complexe, General, ColUpTriang> H(m+1,m+1); H.Fill(zero);

    // s is the vector of residual norm for each inner iteration
    // w is used in the Arnoldi algorithm
    // u is a temporary vector which contains the product A*v(i)
    // r is the residual
    Vector1 w(b), r(b), u(b);
    Vector<Complexe> s(m+1);
    s.Fill(zero); w.Fill(zero); r.Fill(zero); u.Fill(zero);

    for (int i = 0; i < m+1; i++)
      V[i].Fill(zero);

    Vector<Complexe> rotations_sin(m+1);
    rotations_sin.Fill(zero);
    Vector<Titer> rotations_cos(m+1);
    rotations_cos.Fill(Titer(0));

    // we compute residual
    Copy(b, w);
    if (!outer.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), w);
    else
      x.Fill(zero);

    // preconditioning
    M.Solve(A, w, r);
    Titer beta = Norm2(r);

    // we initialize outer
    int success_init = outer.Init(r);
    if (success_init != 0)
      return outer.ErrorCode();

    // the coefficient H(m+1,m)
    Complexe hi_ip1;

    outer.SetNumberIteration(0);
    // Loop until the stopping criteria are reached
    while (! outer.Finished(beta))
      {
	// we normalize V(0) and we init s
	Copy(r, V[0]);
	Mlt(Complexe(Complexe(1)/beta), V[0]);
	s.Fill(zero);
	s(0) = beta;

	int i = 0, k;

	// we initialize the iter iteration
	// m is the maximum number of inner iterations
	Iteration<Titer> inner(outer);
	inner.SetNumberIteration(outer.GetNumberIteration());
	inner.SetMaxNumberIteration(outer.GetNumberIteration()+m);

	do
	  {
	    // product matrix vector u=A*V(i)
	    Mlt(A, V[i], u);

	    // preconditioning
	    M.Solve(A, u, w);

	    // Arnoldi algorithm
	    for (k = 0; k <= i; k++)
	      {
		// h_{k,i} = \bar{v(k)} w
		H.Val(k, i) = DotProdConj(V[k], w);
		Add(-H(k,i), V[k], w);
	      }

	    // we compute h(i+1,i)
	    hi_ip1 = Norm2(w);
	    Copy(w, V[i+1]);

	    // we normalize V(i+1)
	    if (hi_ip1 != zero)
	      Mlt(Complexe(1)/hi_ip1, V[i+1]);

	    // we apply precedent generated rotations
	    // to the last column we computed.
	    for (k = 0; k < i; k++)
	      ApplyRot(H.Val(k,i), H.Val(k+1,i),
		       rotations_cos(k), rotations_sin(k));

	    // we generate a new rotation Omega=[c s;-conj(s) c] in order to
	    // cancel h(i+1,i) and we store this rotation
	    if (hi_ip1 != zero)
	      {
		GenRot(H.Val(i,i), hi_ip1,
		       rotations_cos(i), rotations_sin(i));
	        // After this call we must have hi_ip1=0
		// GenRot must modify the entries H(i,i) hi_ip1
		// we apply the rotation to the right hand side s
		ApplyRot(s(i), s(i+1), rotations_cos(i), rotations_sin(i));
	      }

	    ++inner, ++outer, ++i;

	  } while (! inner.Finished(abs(s(i))));

	// Now we solve the triangular system H
	H.Resize(i, i); s.Resize(i);
	Solve(H, s); s.Resize(m+1);

	// new iterate x = x + sum_0^{i-1} s(k)*V(k)
	for (k = 0; k < i; k++)
	  Add(s(k), V[k], x);

	// we compute the new residual
	Copy(b, w);
	MltAdd(Complexe(-1), A, x, Complexe(1), w);
	M.Solve(A, w, r);

	// residual norm
	beta = Norm2(r);
      }

    return outer.ErrorCode();

  }

} // end namespace

#define SELDON_FILE_ITERATIVE_GMRES_CXX
#endif
