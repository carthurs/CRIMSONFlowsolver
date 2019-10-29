// Copyright (C) 2001-2010 Vivien Mallet
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


#include <cppunit/extensions/HelperMacros.h>

#include "Seldon.hxx"
#include "SeldonSolver.hxx"
using namespace Seldon;


class SparseMatrixTest: public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(SparseMatrixTest);
  CPPUNIT_TEST(test_set_identity);
  CPPUNIT_TEST(test_get_rowcol);
  CPPUNIT_TEST(test_conversion);
  CPPUNIT_TEST(test_permutation);
  CPPUNIT_TEST(test_transposition);
  CPPUNIT_TEST(test_set_rowcol);
  CPPUNIT_TEST_SUITE_END();

protected:
  int m_;
  int n_;
  int Nelement_;
  int Nloop_;

public:
  void setUp()
  {
  }


  void tearDown()
  {
  }


  void test_set_identity()
  {
    m_ = 25;
    n_ = 10;
    set_identity();

    m_ = 10;
    n_ = 25;
    set_identity();

    m_ = 20;
    n_ = 20;
    set_identity();
  }


  void test_get_rowcol()
  {
    srand(time(NULL));

    m_ = 25;
    n_ = 10;
    Nelement_ = 30;
    Nloop_ = 10;
    get_rowcol();

    m_ = 10;
    n_ = 25;
    Nelement_ = 1;
    Nloop_ = 10;
    get_rowcol();

    m_ = 20;
    n_ = 20;
    Nelement_ = 300;
    Nloop_ = 10;
    get_rowcol();
  }


  void test_conversion()
  {
    m_ = 25;
    n_ = 10;
    Nelement_ = 30;
    Nloop_ = 10;
    conversion();

    m_ = 50;
    n_ = 60;
    Nelement_ = 5;
    Nloop_ = 10;
    conversion();

    m_ = 50;
    n_ = 50;
    Nelement_ = 1;
    Nloop_ = 10;
    conversion();

    m_ = 10;
    n_ = 25;
    Nelement_ = 30;
    Nloop_ = 10;
    conversion();

    m_ = 10;
    n_ = 5;
    Nelement_ = 40;
    Nloop_ = 10;
    conversion();

    m_ = 5;
    n_ = 10;
    Nelement_ = 100;
    Nloop_ = 2;
    conversion();
  }


  void test_permutation()
  {
    m_ = 25;
    n_ = 10;
    Nelement_ = 30;
    Nloop_ = 10;
    permutation();

    m_ = 50;
    n_ = 60;
    Nelement_ = 5;
    Nloop_ = 10;
    permutation();

    m_ = 50;
    n_ = 50;
    Nelement_ = 1;
    Nloop_ = 10;
    permutation();

    m_ = 10;
    n_ = 25;
    Nelement_ = 30;
    Nloop_ = 10;
    permutation();

    m_ = 10;
    n_ = 5;
    Nelement_ = 40;
    Nloop_ = 10;
    permutation();

    m_ = 5;
    n_ = 10;
    Nelement_ = 100;
    Nloop_ = 2;
    permutation();
  }


  void test_transposition()
  {
    m_ = 25;
    n_ = 10;
    Nelement_ = 30;
    Nloop_ = 10;
    transposition();

    m_ = 50;
    n_ = 60;
    Nelement_ = 5;
    Nloop_ = 10;
    transposition();

    m_ = 50;
    n_ = 50;
    Nelement_ = 1;
    Nloop_ = 10;
    transposition();

    m_ = 10;
    n_ = 25;
    Nelement_ = 30;
    Nloop_ = 10;
    transposition();

    m_ = 10;
    n_ = 5;
    Nelement_ = 40;
    Nloop_ = 10;
    transposition();

    m_ = 5;
    n_ = 10;
    Nelement_ = 100;
    Nloop_ = 2;
    transposition();
  }


  void test_set_rowcol()
  {
    m_ = 15;
    n_ = 5;
    Nelement_ = 20;
    Nloop_ = 10;
    set_rowcol();

    m_ = 50;
    n_ = 60;
    Nelement_ = 5;
    Nloop_ = 10;
    set_rowcol();

    m_ = 50;
    n_ = 50;
    Nelement_ = 1;
    Nloop_ = 10;
    set_rowcol();

    m_ = 10;
    n_ = 25;
    Nelement_ = 30;
    Nloop_ = 10;
    set_rowcol();

    m_ = 10;
    n_ = 5;
    Nelement_ = 40;
    Nloop_ = 10;
    set_rowcol();
  }


  void set_identity()
  {
    Matrix<double, General, ColSparse> A_col(m_, n_);
    Matrix<double, General, RowSparse> A_row(m_, n_);
    Matrix<double> A_full(m_, n_);

    A_col.SetIdentity();
    A_row.SetIdentity();
    A_full.SetIdentity();

    for (int i = 0; i < m_; i++)
      for (int j = 0; j < n_; j++)
        if (i == j)
        {
          CPPUNIT_ASSERT(A_col(i, j) == double(1));
          CPPUNIT_ASSERT(A_row(i, j) == double(1));
          CPPUNIT_ASSERT(A_full(i, j) == double(1));
        }
        else
        {
          CPPUNIT_ASSERT(A_col(i, j) == double(0));
          CPPUNIT_ASSERT(A_row(i, j) == double(0));
          CPPUNIT_ASSERT(A_full(i, j) == double(0));
        }

    // Testing the function applied to a non-empty matrix.
    A_row.FillRand(m_ + n_);
    A_full.FillRand();
    A_row.SetIdentity();
    A_full.SetIdentity();

    for (int i = 0; i < m_; i++)
      for (int j = 0; j < n_; j++)
        if (i == j)
        {
          CPPUNIT_ASSERT(A_row(i, j) == double(1));
          CPPUNIT_ASSERT(A_full(i, j) == double(1));
        }
        else
        {
          CPPUNIT_ASSERT(A_row(i, j) == double(0));
          CPPUNIT_ASSERT(A_full(i, j) == double(0));
        }
  }


  void get_rowcol()
  {
    int i, j, k;
    for (k = 0; k < Nloop_; k++)
      {
        Matrix<double, General, RowSparse> A(m_, n_);
        A.FillRand(Nelement_);

        /*** To full vectors ***/

        Vector<double> row, column;
        // Rows.
        for (i = 0; i < m_; i++)
          {
            GetRow(A, i, row);
            for (j = 0; j < n_; j++)
              CPPUNIT_ASSERT(A(i, j) == row(j));
          }
        // Columns.
        for (j = 0; j < n_; j++)
          {
            GetCol(A, j, column);
            for (i = 0; i < m_; i++)
              CPPUNIT_ASSERT(A(i, j) == column(i));
          }

        /*** To sparse vectors ***/

        Vector<double, VectSparse> row_sparse, column_sparse;
        // Rows.
        for (i = 0; i < m_; i++)
          {
            GetRow(A, i, row_sparse);
            for (j = 0; j < n_; j++)
              CPPUNIT_ASSERT(A(i, j) == row_sparse(j));
          }
        // Columns.
        for (j = 0; j < n_; j++)
          {
            GetCol(A, j, column_sparse);
            for (i = 0; i < m_; i++)
              CPPUNIT_ASSERT(A(i, j) == column_sparse(i));
          }
      }
  }


  void conversion()
  {
    srand(time(NULL));

    int i, j;
    double value;

    for (int k = 0; k < Nloop_; k++)
      {
        Matrix<double> A_full(m_, n_);
        A_full.Zero();

        Matrix<double, General, ArrayRowSparse> A_array(m_, n_);
        for (int l = 0; l < Nelement_; l++)
          {
            i = rand() % m_;
            j = rand() % n_;
            value = double(rand());
            A_array.AddInteraction(i, j, value);
            A_full(i, j) += value;
          }

        Matrix<double, General, RowSparse> A;

        Copy(A_array, A);

        for (int i = 0; i < m_; i++)
          for (int j = 0; j < n_; j++)
            {
              CPPUNIT_ASSERT(A_full(i, j) == A(i, j));
              CPPUNIT_ASSERT(A_full(i, j) == A_array(i, j));
            }

        Matrix<double, General, ColSparse> A_col;

        Copy(A, A_col);

        for (int i = 0; i < m_; i++)
          for (int j = 0; j < n_; j++)
            CPPUNIT_ASSERT(A_full(i, j) == A_col(i, j));

        Vector<int> row_index, col_index;
        Vector<double> value;
        ConvertMatrix_to_Coordinates(A_array, row_index, col_index, value);

        ConvertMatrix_from_Coordinates(row_index, col_index, value, A);

        for (int i = 0; i < m_; i++)
          for (int j = 0; j < n_; j++)
            CPPUNIT_ASSERT(A_full(i, j) == A(i, j));
      }
  }


  void permutation()
  {
    int i, j;
    double value;

    for (int k = 0; k < Nloop_; k++)
      {
        Matrix<double> A_full(m_, n_);
        A_full.Zero();

        Matrix<double, General, ArrayRowSparse> A_array(m_, n_);
        for (int l = 0; l < Nelement_; l++)
          {
            i = rand() % m_;
            j = rand() % n_;
            value = double(rand());
            A_array.AddInteraction(i, j, value);
            A_full(i, j) += value;
          }

        Matrix<double, General, RowSparse> A;
        Copy(A_array, A);
        Matrix<double, General, ColSparse> A_col;
        Copy(A_array, A_col);

        Vector<int> row_permutation(m_);
        row_permutation.Fill();
        int tmp;
        for (int l = 0; l < m_; l++)
          {
            i = rand() % m_;
            j = rand() % m_;
            tmp = row_permutation(i);
            row_permutation(i) = row_permutation(j);
            row_permutation(j) = tmp;
          }

        Vector<int> col_permutation(n_);
        col_permutation.Fill();
        for (int l = 0; l < n_; l++)
          {
            i = rand() % n_;
            j = rand() % n_;
            tmp = col_permutation(i);
            col_permutation(i) = col_permutation(j);
            col_permutation(j) = tmp;
          }

        ApplyInversePermutation(A, row_permutation, col_permutation);
        ApplyInversePermutation(A_col, row_permutation, col_permutation);

        for (i = 0; i < m_; i++)
          for (j = 0; j < n_; j++)
            {
              CPPUNIT_ASSERT(A_full(i, j) == A(row_permutation(i),
                                               col_permutation(j)));
              CPPUNIT_ASSERT(A_full(i, j) == A_col(row_permutation(i),
                                                   col_permutation(j)));
            }
      }
  }


  void transposition()
  {
    srand(time(NULL));

    int i, j;
    double value;

    for (int k = 0; k < Nloop_; k++)
      {
	{
	  Matrix<double, General, ArrayRowSparse> A_array(m_, n_),
            A_array_t(m_, n_);
	  for (int l = 0; l < Nelement_; l++)
	    {
	      i = rand() % m_;
	      j = rand() % n_;
	      value = double(rand());
	      A_array.AddInteraction(i, j, value);
	    }

	  Copy(A_array, A_array_t);
	  Transpose(A_array_t);

	  Matrix<double, General, RowSparse, MallocAlloc<double> > A(m_, n_);
	  Copy(A_array, A);

	  Transpose(A);

	  for (int i = 0; i < m_; i++)
	    for (int j = 0; j < n_; j++)
	      CPPUNIT_ASSERT(A_array_t(j, i) == A(j, i));

	  Transpose(A);

	  for (int i = 0; i < m_; i++)
	    for (int j = 0; j < n_; j++)
	      CPPUNIT_ASSERT(A_array(i, j) == A(i, j));
	}

	{
	  Matrix<double, General, ArrayRowSparse> A_array(m_, n_),
            A_array_t(m_, n_);
	  for (int l = 0; l < Nelement_; l++)
	    {
	      i = rand() % m_;
	      j = rand() % n_;
	      value = double(rand());
	      A_array.AddInteraction(i, j, value);
	    }

	  Copy(A_array, A_array_t);
	  Transpose(A_array_t);

	  Matrix<double, General, RowSparse, CallocAlloc<double> > A(m_, n_);
	  Copy(A_array, A);

	  Transpose(A);

	  for (int i = 0; i < m_; i++)
	    for (int j = 0; j < n_; j++)
	      CPPUNIT_ASSERT(A_array_t(j, i) == A(j, i));

	  Transpose(A);

	  for (int i = 0; i < m_; i++)
	    for (int j = 0; j < n_; j++)
	      CPPUNIT_ASSERT(A_array(i, j) == A(i, j));
	}

	{
	  Matrix<double, General, ArrayRowSparse> A_array(m_, n_),
            A_array_t(m_, n_);
	  for (int l = 0; l < Nelement_; l++)
	    {
	      i = rand() % m_;
	      j = rand() % n_;
	      value = double(rand());
	      A_array.AddInteraction(i, j, value);
	    }

	  Copy(A_array, A_array_t);
	  Transpose(A_array_t);

	  Matrix<double, General, RowSparse, MallocObject<double> > A(m_, n_);
	  Copy(A_array, A);

	  Transpose(A);

	  for (int i = 0; i < m_; i++)
	    for (int j = 0; j < n_; j++)
	      CPPUNIT_ASSERT(A_array_t(j, i) == A(j, i));

	  Transpose(A);

	  for (int i = 0; i < m_; i++)
	    for (int j = 0; j < n_; j++)
	      CPPUNIT_ASSERT(A_array(i, j) == A(i, j));
	}

	{
	  Matrix<double, General, ArrayRowSparse> A_array(m_, n_),
            A_array_t(m_, n_);
	  for (int l = 0; l < Nelement_; l++)
	    {
	      i = rand() % m_;
	      j = rand() % n_;
	      value = double(rand());
	      A_array.AddInteraction(i, j, value);
	    }

	  Copy(A_array, A_array_t);
	  Transpose(A_array_t);

	  Matrix<double, General, RowSparse, NewAlloc<double> > A(m_, n_);
	  Copy(A_array, A);

	  Transpose(A);

	  for (int i = 0; i < m_; i++)
	    for (int j = 0; j < n_; j++)
	      CPPUNIT_ASSERT(A_array_t(j, i) == A(j, i));

	  Transpose(A);

	  for (int i = 0; i < m_; i++)
	    for (int j = 0; j < n_; j++)
	      CPPUNIT_ASSERT(A_array(i, j) == A(i, j));
	}
      }
  }


  void set_rowcol()
  {
    srand(time(NULL));

    int i, j, k, Nrow_sparse, Ncol_sparse;

    for (k = 0; k < Nloop_; k++)
      {
	{
	  Matrix<double, General, RowSparse, MallocAlloc<double> > A(m_, n_);
	  A.FillRand(Nelement_);
	  for (i = 0; i < m_; i++)
	    {
	      Nrow_sparse =  rand() % n_ + 1;
	      Vector<double, VectSparse> row_sparse(Nrow_sparse);
	      for (int l = 0; l < Nrow_sparse; l++)
		row_sparse.Index(l) = rand() % n_;
	      row_sparse.FillRand();
	      row_sparse.Assemble();

	      SetRow(row_sparse, i, A);

	      for (j = 0; j < n_; j++)
		CPPUNIT_ASSERT(A(i, j) == row_sparse(j));
	    }
	}
	{
	  Matrix<double, General, RowSparse, MallocAlloc<double> > A(m_, n_);
	  A.FillRand(Nelement_);
	  for (j = 0; j < n_; j++)
	    {
	      Ncol_sparse =  rand() % m_ + 1;
	      Vector<double, VectSparse> col_sparse(Ncol_sparse);
	      for (int l = 0; l < Ncol_sparse; l++)
		col_sparse.Index(l) = rand() % m_;
	      col_sparse.FillRand();
	      col_sparse.Assemble();

	      SetCol(col_sparse, j, A);

	      for (i = 0; i < m_; i++)
		CPPUNIT_ASSERT(A(i, j) == col_sparse(i));
	    }
	}
      }
  }

};
