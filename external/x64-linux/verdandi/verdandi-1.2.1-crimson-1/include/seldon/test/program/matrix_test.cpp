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


#define SELDON_DEBUG_LEVEL_4
#define SELDON_DEFAULT_ALLOCATOR NewAlloc

#include "Seldon.hxx"

using namespace Seldon;

// testing class Matrix without Blas

template<class T, class Prop, class Storage, class Allocator>
void TestGeneralMatrix(Matrix<T, Prop, Storage, Allocator>& mat_test)
{
  // Test of General dense matrices
  Matrix<T, Prop, Storage, Allocator> A(3,2), B;
  B.Reallocate(4,4);

  // A and B are not initialized
  // you can initialize them with Zero and SetIdentity for example
  A.Zero();
  B.SetIdentity();

  cout << "Number of rows in A : " << A.GetM() << endl;
  cout << "Number of columns in A : " << A.GetN() << endl;
  cout << "Number of elements in A : " << A.GetSize() << endl;
  cout << " A = " << endl << A << endl; // should return zero matrix
  cout << " B = " << endl << B << endl; // should return identity matrix

  // Reallocate doesn't keep previous elements
  B.Reallocate(3,3);
  // so you can initialize B with Fill for example
  B.FillRand();
  // Resize keeps previous elements
  A.Fill();
  cout << " A = " << A << endl;
  A.Resize(2, 4);
  // but new elements need to be initialized
  A(0, 2) = -2.5; A(0, 3) = 1.0;
  A(1, 2) = 3.2; A(1, 3) = -0.5;
  cout << " A = " << A << endl;

  // you can clear a matrix
  A.Clear();
  cout << "Number of elements in A : " << A.GetDataSize() << endl;
  A.Reallocate(3, 3);
  A.Fill(2.2);
  Matrix<T, Prop, Storage, Allocator> C(A);
  cout << "C = " << C << endl;
  // operators = and *= are available
  C = A;
  C *= -2.5;
  cout << "C = " << C << endl;
  // = can be used as an equivalent of method Fill
  C = 1.5;
  cout << "C = " << C << endl;

  // you can use functions GetData, GetIndex, SetData and Nullify
  // for low level manipulations
  Matrix<T, Prop, Storage, Allocator> num(4,2), num2;
  num.Fill();
  T* data = num.GetData();
  // you set num2 with pointer data
  num2.SetData(4, 2, data);
  // in order to avoid that the destructor of num and num2
  // release both the same memory, you can use nullify
  num.Nullify();  // num is now empty

  // you can write and read matrices in files
  A.Fill();
  A.Write("mat_binary.dat");
  A.WriteText("mat_ascii.dat");
  B.Clear();
  B.Read("mat_binary.dat");
  cout << " B = " << B << endl;
  B.ReadText("mat_ascii.dat");
  cout << " B = " << B << endl;

  // a row of A is extracted
  Vector<T> rowA;
  GetRow(A, 1, rowA);
  cout << "row 1 of A = " << rowA << endl;
  // you can change a row of A
  rowA.Fill(2.3);
  SetRow(rowA, 2, A);
  cout << "A = " << A << endl;

  // same stuff for columns
  Vector<T> colA;
  GetCol(A, 0, colA);
  cout << "column 0 of A = " << colA << endl;
  colA.Fill(-1.1);
  SetCol(colA, 1, A);
  cout << "A = " << A << endl;
}


template<class T, class Prop, class Storage, class Allocator>
void TestSymmetricMatrix(Matrix<T, Prop, Storage, Allocator>& mat_test)
{
  // Test of General dense matrices
  Matrix<T, Prop, Storage, Allocator> A(3,3), B;
  B.Reallocate(4,4);

  // A and B are not initialized
  // you can initialize them with Zero and SetIdentity for example
  A.Zero();
  B.SetIdentity();

  cout << "Number of rows in A : " << A.GetM() << endl;
  cout << "Number of columns in A : " << A.GetN() << endl;
  cout << "Number of elements in A : " << A.GetSize() << endl;

  cout << " B = " << endl << B << endl; // should return identity matrix

  // Reallocate doesn't keep previous elements
  B.Reallocate(3,3);
  // so you can initialize B with Fill for example
  B.FillRand();
  // Resize keeps previous elements
  A.Fill(); A.Resize(4, 4);
  // but new elements need to be initialized
  // symmetric matrices -> only the last column need to be changed
  A.Get(0, 3) = -2.5; A.Get(1, 3) = 1.0;
  A.Get(2, 3) = 3.0; A.Get(3, 3) = -0.5;
  cout << " A = " << A << endl;

  // you can clear a matrix
  A.Clear();
  cout << "Number of elements in A : " << A.GetDataSize() << endl;
  A.Reallocate(3, 3);
  A.Fill(-4.1);
  Matrix<T, Prop, Storage, Allocator> C(A);
  // operators = and *= are available
  C = A;
  C *= -2.5;
  cout << "C = " << C << endl;
  // = can be used as an equivalent of method Fill
  C = 1.5;
  cout << "C = " << C << endl;

  // you can use functions GetData, GetIndex, SetData and Nullify
  // for low level manipulations
  Matrix<T, Prop, Storage, Allocator> num(4,4), num2;
  num.Fill();
  T* data = num.GetData();
  // you set num2 with pointer data
  num2.SetData(4, 4, data);
  // in order to avoid that the destructor of num and num2
  // release both the same memory, you can use nullify
  num.Nullify();  // num is now empty

  // you can write and read matrices in files
  A.Fill();
  A.Write("mat_binary.dat");
  A.WriteText("mat_ascii.dat");
  B.Clear();
  B.Read("mat_binary.dat");
  cout << " B = " << B << endl;
  B.ReadText("mat_ascii.dat");
  cout << " B = " << B << endl;

}


template<class T, class Prop, class Storage, class Allocator>
void TestHermitianMatrix(Matrix<T, Prop, Storage, Allocator>& mat_test)
{
  // Test of General dense matrices
  Matrix<T, Prop, Storage, Allocator> A(3,3), B;
  B.Reallocate(4,4);

  // A and B are not initialized
  // you can initialize them with Zero and SetIdentity for example
  A.Zero();
  B.SetIdentity();

  cout << "Number of rows in A : " << A.GetM() << endl;
  cout << "Number of columns in A : " << A.GetN() << endl;
  cout << "Number of elements in A : " << A.GetSize() << endl;

  cout << " B = " << endl << B << endl; // should return identity matrix

  // Reallocate doesn't keep previous elements
  B.Reallocate(3,3);
  // so you can initialize B with Fill for example
  B.FillRand();
  // Resize keeps previous elements
  A.Fill(); A.Resize(4, 4);
  // but new elements need to be initialized
  // hermitian matrices -> only the last column need to be changed
  // be careful, we don't check if the diagonal is real
  A.Get(0, 3) = T(-2.5, 1.0); A.Get(1, 3) = T(1.0, -1.5);
  A.Get(2, 3) = T(3.0, 2.0); A.Get(3, 3) = T(-0.5, 0.0);
  cout << " A = " << A << endl;

  // you can clear a matrix
  A.Clear();
  cout << "Number of elements in A : " << A.GetDataSize() << endl;
  A.Reallocate(3, 3);
  A.Fill(-1.9);
  Matrix<T, Prop, Storage, Allocator> C(A);
  cout << " C = " << C << endl;
  // operators = and *= are available
  C = A;
  C *= T(-2.5, 1.8);
  cout << "C = " << C << endl;
  // = can be used as an equivalent of method Fill
  C = 1.5;
  cout << "C = " << C << endl;

  // you can use functions GetData, GetIndex, SetData and Nullify
  // for low level manipulations
  Matrix<T, Prop, Storage, Allocator> num(4,4), num2;
  num.Fill();
  T* data = num.GetData();
  // you set num2 with pointer data
  num2.SetData(4, 4, data);
  // in order to avoid that the destructor of num and num2
  // release both the same memory, you can use nullify
  num.Nullify();  // num is now empty

  // you can write and read matrices in files
  A.Fill();
  A.Write("mat_binary.dat");
  A.WriteText("mat_ascii.dat");
  B.Clear();
  B.Read("mat_binary.dat");
  cout << " B = " << B << endl;
  B.ReadText("mat_ascii.dat");
  cout << " B = " << B << endl;

}


template<class T, class Prop, class Storage, class Allocator>
void TestUpperTriangularMatrix(Matrix<T, Prop, Storage, Allocator>& mat_test)
{
  // Test of General dense matrices
  Matrix<T, Prop, Storage, Allocator> A(3,3), B;
  B.Reallocate(4,4);

  // A and B are not initialized
  // you can initialize them with Zero and SetIdentity for example
  A.Zero();
  B.SetIdentity();

  cout << "Number of rows in A : " << A.GetM() << endl;
  cout << "Number of columns in A : " << A.GetN() << endl;
  cout << "Number of elements in A : " << A.GetSize() << endl;

  cout << " B = " << endl << B << endl; // should return identity matrix

  // Reallocate doesn't keep previous elements
  B.Reallocate(3,3);
  // so you can initialize B with Fill for example
  B.FillRand();
  // Resize keeps previous elements
  A.Fill(); A.Resize(4, 4);
  // but new elements need to be initialized
  // upper triangular matrices -> only the last column need to be changed
  A.Get(0, 3) = -2.5; A.Get(1, 3) = 1.0;
  A.Get(2, 3) = 3.2; A.Get(3, 3) = -0.5;
  cout << " A = " << A << endl;

  // you can clear a matrix
  A.Clear();
  cout << "Number of elements in A : " << A.GetDataSize() << endl;
  A.Reallocate(3, 3);
  A.Fill(-1.9);
  Matrix<T, Prop, Storage, Allocator> C(A);
  // operators = and *= are available
  C = A;
  C *= -2.5;
  cout << "C = " << C << endl;
  // = can be used as an equivalent of method Fill
  C = 1.5;
  cout << "C = " << C << endl;

  // you can use functions GetData, GetIndex, SetData and Nullify
  // for low level manipulations
  Matrix<T, Prop, Storage, Allocator> num(4,4), num2;
  num.Fill();
  T* data = num.GetData();
  // you set num2 with pointer data
  num2.SetData(4, 4, data);
  // in order to avoid that the destructor of num and num2
  // release both the same memory, you can use nullify
  num.Nullify();  // num is now empty

  // you can write and read matrices in files
  A.Fill();
  A.Write("mat_binary.dat");
  A.WriteText("mat_ascii.dat");
  B.Clear();
  B.Read("mat_binary.dat");
  cout << " B = " << B << endl;
  B.ReadText("mat_ascii.dat");
  cout << " B = " << B << endl;

}


template<class T, class Prop, class Storage, class Allocator>
void TestLowerTriangularMatrix(Matrix<T, Prop, Storage, Allocator>& mat_test)
{
  // Test of General dense matrices
  Matrix<T, Prop, Storage, Allocator> A(3,3), B;
  B.Reallocate(4,4);

  // A and B are not initialized
  // you can initialize them with Zero and SetIdentity for example
  A.Zero();
  B.SetIdentity();

  cout << "Number of rows in A : " << A.GetM() << endl;
  cout << "Number of columns in A : " << A.GetN() << endl;
  cout << "Number of elements in A : " << A.GetSize() << endl;

  cout << " B = " << endl << B << endl; // should return identity matrix

  // Reallocate doesn't keep previous elements
  B.Reallocate(3,3);
  // so you can initialize B with Fill for example
  B.FillRand();
  // Resize keeps previous elements
  A.Fill(); A.Resize(4, 4);
  // but new elements need to be initialized
  // lower triangular matrices -> only the last row need to be changed
  A.Val(3, 0) = -2.5; A.Val(3, 1) = 1.0;
  A.Val(3, 2) = 3.2; A.Val(3, 3) = -0.5;
  cout << " A = " << A << endl;

  // you can clear a matrix
  A.Clear();
  cout << "Number of elements in A : " << A.GetDataSize() << endl;
  A.Reallocate(3, 3);
  A.Fill(-1.9);
  Matrix<T, Prop, Storage, Allocator> C(A);
  // operators = and *= are available
  C = A;
  C *= -2.5;
  cout << "C = " << C << endl;
  // = can be used as an equivalent of method Fill
  C = 1.5;
  cout << "C = " << C << endl;

  // you can use functions GetData, GetIndex, SetData and Nullify
  // for low level manipulations
  Matrix<T, Prop, Storage, Allocator> num(4,4), num2;
  num.Fill();
  T* data = num.GetData();
  // you set num2 with pointer data
  num2.SetData(4, 4, data);
  // in order to avoid that the destructor of num and num2
  // release both the same memory, you can use nullify
  num.Nullify();  // num is now empty

  // you can write and read matrices in files
  A.Fill();
  A.Write("mat_binary.dat");
  A.WriteText("mat_ascii.dat");
  B.Clear();
  B.Read("mat_binary.dat");
  cout << " B = " << B << endl;
  B.ReadText("mat_ascii.dat");
  cout << " B = " << B << endl;

}


int main()
{
  cout << "Seldon: compilation test of class Matrix without Blas" << endl;

  Matrix<double, General, RowMajor> A1;
  TestGeneralMatrix(A1);

  Matrix<double, General, ColMajor> A2;
  TestGeneralMatrix(A2);

  //Matrix<double, Symmetric, RowSym> A3;
  //TestSymmetricMatrix(A3);

  //Matrix<double, Symmetric, ColSym> A4;
  //TestSymmetricMatrix(A4);

  Matrix<complex<double>, Symmetric, RowSymPacked> A5;
  TestSymmetricMatrix(A5);

  Matrix<complex<double>, Symmetric, ColSymPacked> A6;
  TestSymmetricMatrix(A6);

  //Matrix<complex<double>, General, RowHerm> A7
  //TestHermitianMatrix(A7);

  //Matrix<complex<double>, General, ColHerm> A8;
  //TestHermitianMatrix(A8);

  Matrix<complex<double>, General, RowHermPacked> A9;
  TestHermitianMatrix(A9);

  Matrix<complex<double>, General, ColHermPacked> A10;
  TestHermitianMatrix(A10);

  //Matrix<double, General, RowLoTriang> A11;
  //TestLowerTriangularMatrix(A11);

  //Matrix<double, General, ColLoTriang> A12;
  //TestLowerTriangularMatrix(A12);

  Matrix<double, General, RowLoTriangPacked> A13;
  TestLowerTriangularMatrix(A13);

  Matrix<double, General, ColLoTriangPacked> A14;
  TestLowerTriangularMatrix(A14);

  //Matrix<double, General, RowUpTriang> A15;
  //TestLowerTriangularMatrix(A15);

  //Matrix<double, General, ColUpTriang> A16;
  //TestLowerTriangularMatrix(A16);

  Matrix<double, General, RowUpTriangPacked> A17;
  TestUpperTriangularMatrix(A17);

  Matrix<double, General, ColUpTriangPacked> A18;
  TestUpperTriangularMatrix(A18);

  return 0;
}
