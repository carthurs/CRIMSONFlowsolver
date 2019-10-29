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


#define SELDON_DEFAULT_ALLOCATOR NewAlloc

#include "Seldon.hxx"

using namespace Seldon;

// testing class Vector without Blas

int main()
{
  cout << "Seldon: compilation test of class Vector without Blas" << endl;

  Vector<double> U(3), V;
  V.Reallocate(3);
  // U and V are not initialized
  // you can initialize them with Zero, Fill or FillRand
  U.Zero();
  cout << " U = " << U << endl;
  V.Fill();
  cout << "Number of elements in U : " << U.GetM() << endl;
  cout << "Number of elements in V : " << V.GetLength() << endl;
  cout << " V : " << V << endl; // should return [0 1 2]

  // Reallocate doesn't keep previous elements
  V.Reallocate(5);
  cout << "Number of elements in V : " << V.GetSize() << endl;
  // so you can initialize V
  V.FillRand();
  // Resize keeps previous elements
  U.Fill(1.5); U.Resize(5);
  // but new elements need to be initialized
  U(3) = 1.2;
  U(4) = -1.0;
  cout << "Vector U : " << U << endl;

  // you can erase a vector
  V.Clear();
  cout << "Number of elements in V : " << V.GetDataSize() << endl;
  // use a copy constructor
  V.Reallocate(2);
  V.Fill();
  Vector<double> W = V;
  cout << "Vector W : " << W << endl;
  // operators = and *= are available
  W = U;
  W *= -2.5;
  cout << "Vector W : " << W << endl;
  // = can be used as an equivalent of method Fill
  W = 1.5;
  cout << "Vector W : " << W << endl;

  // you can use functions GetData, GetIndex, SetData and Nullify
  // for low level manipulations
  IVect num(4), num2;
  num.Fill();
  int* data = num.GetData();
  // you set num2 with pointer data
  num2.SetData(4, data);
  cout << " num2 = " << num2 << endl;
  // in order to avoid that the destructor of num and num2
  // release both the same memory, you can use nullify
  num.Nullify();  // num is now empty

  // you can insert elements to the end with PushBack
  Vector<string> A, B;
  A.PushBack(string("flower"));
  B.Reallocate(2);
  B(0) = "animal"; B(1) = "human";
  A.PushBack(B);
  cout<<" A = "; A.Print(); cout << endl;

  // GetNormInfIndex returns the index
  // where the absolute highest value is reached
  num.Reallocate(3);
  num(0) = 0; num(1) = -5; num(2) = 3;
  int imax = num.GetNormInfIndex();
  cout << "Index where |num| reaches its maximum : " << imax <<endl;

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
    Vector<double> X(4), Y(4);

    X.Fill();
    Y.Fill(1.4);

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
    Vector<complex<double> > X(4), Y(4);

    X.Fill();
    X(2) = complex<double>(1, -1);
    Y.Fill(complex<double>(1.4, 0.2) );

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

  // you can sort a vector
  Vector<int> vec(5), vec2;
  vec(0) = 1; vec(1) = 2; vec(2) = 0; vec(3) = -3; vec(4) = 1;
  vec2 = vec;
  Sort(vec2);
  // the result should be -3 0 1 1 2
  cout << "Sorted vector " << vec2 << endl;
  // you can sort only a part of the vector
  vec2 = vec;
  Sort(1, 3, vec2);
  // the result should be 1 -3 0 2 1
  cout << "Vector partially sorted " << vec2 << endl;
  // you can sort a vector, and the sorting operation affects other arrays
  // this is very useful to retrieve the permutation vector
  IVect permutation(5); permutation.Fill();
  vec2 = vec;
  // permutation should become 3 2 4 0 1  (or 3 2 0 4 1)
  Sort(vec2, permutation);
  cout << "Permutation vector " << permutation << endl;
  // you can sort and remove duplicate elements by calling Assemble
  // be careful, because the resulting vector is not resized
  vec2 = vec;
  int nb = vec2.GetM();
  Assemble(nb, vec2);
  // nb should be equal to 4 whereas vec2.GetM() is equal to 5
  // you can resize the vector, vec2 should then be equal to -3 0 1 2
  vec2.Resize(nb);
  cout << "Sorted vector without duplicate elements " << vec2 << endl;
  // you can add values of a second vector with function
  Vector<double> values(5);
  values(0) = 0.3; values(1) = -0.5;
  values(2) = 0.7; values(3) = 0.4; values(4) = -0.1;
  vec2 = vec; nb = 5;
  Assemble(nb, vec2, values);
  vec2.Resize(nb); values.Resize(nb);
  // values should be equal to 0.4 0.7 0.2 -0.5
  // the third element 0.2 comes from the operation 0.3 - 0.1
  cout << "values " << values << endl;
  // you can remove duplicate of a first vector and affect a second one
  permutation.Fill();
  vec2 = vec; nb = 5;
  RemoveDuplicate(nb, vec2, permutation);
  // again the vectors are not resized
  vec2.Resize(nb); permutation.Resize(nb);
  // permutation should become 3 2 0 1  (or 3 2 4 1)
  cout << "Permutation vector " << permutation << endl;

  return 0;

}
