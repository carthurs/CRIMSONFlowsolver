// Copyright (C) 2010 Marc Duruflé
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


#include "Seldon.hxx"
#include "SeldonSolver.hxx"

using namespace Seldon;

template<class T>
bool CheckVector(Vector<T>& x)
{
  bool test = true;
  for (int i = 0; i < x.GetM(); i++)
    if (abs(x(i) - i)/x.GetM() > 1e-15)
      test = false;

  return test;
}


int main()
{
  TRY;

  cout.precision(16);

  // symmetric matrix is read in a file
  Matrix<double, Symmetric, ArrayRowSymSparse > A;
  A.ReadText("matrix/MhSparse.dat");

  // creation of a right hand side b = A*[0;1;...;n-1]
  Vector<double> x(A.GetM()), b(A.GetM()), y(A.GetM());
  x.Fill();
  Mlt(A, x, b);

  // Cholesky factorization using Seldon function
  GetCholesky(A);

  x = b;
  SolveCholesky(SeldonNoTrans, A, x);
  SolveCholesky(SeldonTrans, A, x);

  bool test_solve_seldon = CheckVector(x);

  y.Fill();
  SolveCholesky(SeldonNoTrans, A, y);
  SolveCholesky(SeldonTrans, A, y);

  x = y;
  MltCholesky(SeldonTrans, A, x);
  MltCholesky(SeldonNoTrans, A, x);

  bool test_mlt_seldon = CheckVector(x);

#ifdef SELDON_WITH_CHOLMOD
  // Cholesky factorization using Cholmod.
  MatrixCholmod mat_chol;

  A.ReadText("matrix/MhSparse.dat");

  GetCholesky(A, mat_chol);

  x = b;
  SolveCholesky(SeldonNoTrans, mat_chol, x);
  SolveCholesky(SeldonTrans, mat_chol, x);

  bool test_solve_cholmod = CheckVector(x);
#endif

  bool all_test = true;
  if (!test_solve_seldon)
    {
      all_test = false;
      cout << "SolveCholesky provided by Seldon incorrect" << endl;
    }

  if (!test_mlt_seldon)
    {
      all_test = false;
      cout << "MltCholesky provided by Seldon incorrect" << endl;
    }

#ifdef SELDON_WITH_CHOLMOD
  if (!test_solve_cholmod)
    {
      all_test = false;
      cout << "SolveCholesky provided by Cholmod incorrect" << endl;
    }
#endif

  A.ReadText("matrix/MatFente.dat");

  // creation of a right hand side b = A*[0;1;...;n-1]
  x.Reallocate(A.GetM());
  b.Reallocate(A.GetM());
  y.Reallocate(A.GetM());
  x.FillRand();
  Mlt(1e-9, x);
  Vector<double> xsol = x;
  Mlt(A, x, b);

  // Cholesky factorization using Seldon function.
  SparseCholeskySolver<double> solver;
  solver.ShowMessages();
  solver.SelectOrdering(SparseMatrixOrdering::METIS);
  solver.SelectDirectSolver(solver.SELDON_SOLVER);
  // solver.SelectDirectSolver(solver.CHOLMOD);
  solver.Factorize(A);

  x = b;
  solver.Solve(SeldonNoTrans, x);
  solver.Solve(SeldonTrans, x);

  for (int i = 0; i < x.GetM(); i++)
    if (abs(x(i) - xsol(i)) > 1e-12)
      {
        cout << "Solver failed." << endl;
        abort();
      }

  solver.Mlt(SeldonTrans, x);
  solver.Mlt(SeldonNoTrans, x);

  for (int i = 0; i < x.GetM(); i++)
    if (abs(x(i) - b(i)) > 1e-12)
      {
        cout << "Mlt failed." << endl;
        abort();
      }


  if (all_test)
    cout << "All tests passed successfully" << endl;
  else
    return -1;

  END;

  return 0;
}
