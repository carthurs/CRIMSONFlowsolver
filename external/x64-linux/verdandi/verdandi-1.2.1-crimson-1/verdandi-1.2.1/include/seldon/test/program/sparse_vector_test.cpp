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

// testing class Vector<T, VectSparse>

int main()
{
  cout << "Seldon: compilation test of sparse vectors" << endl;

  Vector<double, VectSparse> U(3), V;
  V.Reallocate(2); V.Reallocate(3);
  // U and V are not initialized
  // you can initialize them with Zero, Fill or FillRand
  U.Zero();
  V.Fill();
  // but this is not sufficient, since you have to specify row numbers as well
  // you can do that with method Index
  // usually row numbers should be given sorted
  U.Index(0) = 0;
  U.Index(1) = 4;
  U.Index(2) = 7;
  // but you can give them unsorted
  V.Index(0) = 1;
  V.Index(1) = 5;
  V.Index(2) = 3;
  // if row numbers are unsorted, you have to call Assemble
  V.Assemble();
  cout << "Number of elements in U : " << U.GetM() << endl;
  cout << "Number of elements in V : " << V.GetLength() << endl;
  cout << " U : " << U << endl;
  cout << " V : " << V << endl;

  // Reallocate doesn't keep previous elements
  V.Reallocate(5);
  cout << "Number of elements in V : " << V.GetSize() << endl;
  // so you can initialize V
  V.FillRand();
  V.Index(0) = 1;
  V.Index(1) = 2;
  V.Index(2) = 6;
  V.Index(3) = 9;
  V.Index(4) = 11;
  // Resize keeps previous elements
  U.Fill(1.5); U.Resize(5);
  // but new elements need to be initialized
  U.Index(3) = 3; U.Value(3) = 1.2;
  U.Index(4) = 6; U.Value(4) = -1.0;
  // again if new numbers are not sorted, you have to call Assemble
  U.Assemble();
  cout << "Vector U : " << U << endl;

  // you can erase a vector
  V.Clear();
  cout << "Number of elements in V : " << V.GetDataSize() << endl;
  // use a copy constructor
  V.Reallocate(2);
  V.Fill();
  V.Index(0) = 5; V.Index(1) = 10;
  Vector<double, VectSparse> W = V;
  cout << "Vector W : " << W << endl;
  // operators = and *= are available
  W = U;
  W *= -2.5;
  cout << "Vector W : " << W << endl;
  // = can be used as an equivalent of method Fill
  W = 1.5;
  cout << "Vector W : " << W << endl;

  // you can use functions GetData, SetData and Nullify
  // for low level manipulations
  Vector<int> num(4); Vector<double> values(4);
  num(0) = 1; values(0) = 1.1;
  num(1) = 3; values(1) = -0.3;
  num(2) = 5; values(2) = 0.2;
  num(3) = 9; values(3) = 0.7;
  // you set row numbers and values of U
  U.SetData(values, num);
  cout<<" U = " << U << endl;
  // values and num are empty after the call to SetData
  // you can retrieve row numbers and values of W with GetData and GetIndex
  num.SetData(W.GetM(), W.GetIndex());
  values.SetData(W.GetM(), W.GetData());
  // when you modify num and values, W is changed
  num(3) = 5; values(3) = 2.2;
  cout << "Vector W = " << W << endl;
  // in order to avoid that the destructor of W, num and values
  // release both the same memory, you can use nullify
  num.Nullify();  // num is now empty
  values.Nullify(); // values is now empty

  // you can insert elements by using AddInteraction
  Vector<double, VectSparse> A, B;
  A.AddInteraction(4, 1.5);
  cout << " A " << A << endl;
  // or AddInteractionRow for several values
  num.Reallocate(3); values.Reallocate(3);
  num(0) = 7; values(0) = -1.0;
  num(1) = 2; values(1) = 0.2;
  num(2) = 4; values(2) = 0.5;
  // you dont need to sort those values
  // if vector V already contains a value
  // it is added in the same position
  A.AddInteractionRow(num.GetM(), num, values);
  cout<<" A = "; A.Print(); cout << endl;
  // you can also use operator () to display the value
  // of A at a given row number or modify it
  cout << "A(3) = " << A(3) << endl;
  // Use Get if you want to modify the value.
  A.Get(3) = 0.5;
  cout << "A(3) = " << A(3) << endl;

  // GetNormInfIndex returns the index
  // where the absolute highest value is reached
  int imax = A.GetNormInfIndex();
  cout << "Index where |A| reaches its maximum : " << imax <<endl;
  // if you want the row number and absolute highest value
  // don't forget to use Index and Value
  cout << "Maximum " << A.Value(imax)
       << " is reached for row " << A.Index(imax) << endl;

  // you can write and read vectors in files
  U.Write("vec_binary.dat");
  U.WriteText("vec_ascii.dat");
  V.Clear();
  V.Read("vec_binary.dat");
  cout << " V = " << V << endl;
  W.ReadText("vec_ascii.dat");
  cout << " W = " << W << endl;

  // you can use functions Mlt, Add, Copy, Swap, DotProd, Norm1, Norm2
  {
    // for real numbers
    Vector<double, VectSparse> X(4), Y(4);

    X.Index(0) = 2; X.Value(0) = -0.4;
    X.Index(1) = 5; X.Value(1) = 1.6;
    X.Index(2) = 11; X.Value(2) = 0.8;
    X.Index(3) = 12; X.Value(3) = -1.1;

    Y.Index(0) = 0; Y.Value(0) = 1.7;
    Y.Index(1) = 5; Y.Value(1) = -2.2;
    Y.Index(2) = 7; Y.Value(2) = 0.6;
    Y.Index(3) = 12; Y.Value(3) = 0.9;

    Mlt(1.3, X);
    Add(2.2, X, Y);
    cout << " Y = " << Y << endl;

    Swap(X, Y);
    cout << " X = " << X << endl;
    cout << " Y = " << Y << endl;

    Copy(X, Y);
    cout << " Y = " << Y << endl;

    Y.Fill(-1.6);
    double scal = DotProd(X, Y);
    cout << "Scalar product " << scal << endl;

    scal = Norm2(X);
    cout << "2-Norm of X " << scal << endl;

    scal = Norm1(X);
    cout << "1-Norm of X " << scal << endl;
  }

  {
    // for complex numbers
    Vector<complex<double> , VectSparse> X(4), Y(4);

    X.Index(0) = 2; X.Value(0) = complex<double>(-0.4, 0.4);
    X.Index(1) = 5; X.Value(1) = complex<double>(1.6, -0.2);
    X.Index(2) = 11; X.Value(2) = complex<double>(0.8, 0.3);
    X.Index(3) = 12; X.Value(3) = complex<double>(-1.1, 0);

    Y.Index(0) = 0; Y.Value(0) = complex<double>(1.7, -0.1);
    Y.Index(1) = 5; Y.Value(1) = complex<double>(-2.2, 2.8);
    Y.Index(2) = 7; Y.Value(2) = complex<double>(0.6, -1.4);
    Y.Index(3) = 12; Y.Value(3) = complex<double>(0.9, 3.5);

    Mlt(complex<double>(1.3, 0.5), X);
    Add(complex<double>(2.2, -0.7), X, Y);
    cout << " Y = " << Y << endl;

    Swap(X, Y);
    cout << " X = " << X << endl;
    cout << " Y = " << Y << endl;

    Copy(X, Y);
    cout << " Y = " << Y << endl;

    Y.Fill(complex<double>(-1.6, 0.3) );
    complex<double> scal = DotProd(X, Y);
    cout << "Scalar product " << scal << endl;

    double scalb = Norm2(X);
    cout << "2-Norm of X " << scalb << endl;

    scalb = Norm1(X);
    cout << "1-Norm of X " << scalb << endl;

    // vector can be conjugated
    Conjugate(X);
    cout << " X = " << X << endl;
  }

  return 0;

}
