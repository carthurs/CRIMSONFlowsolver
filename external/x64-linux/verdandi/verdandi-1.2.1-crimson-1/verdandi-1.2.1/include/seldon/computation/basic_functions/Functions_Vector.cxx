// Copyright (C) 2003-2009 Marc Duruflé
// Copyright (C) 2001-2009 Vivien Mallet
// Copyright (C) 2010 INRIA
// Author(s): Marc Fragu
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


#ifndef SELDON_FILE_FUNCTIONS_VECTOR_CXX
#define SELDON_FILE_FUNCTIONS_VECTOR_CXX


#include "Functions_Vector.hxx"


/*
  Functions defined in this file:

  alpha X -> X
  Mlt(alpha, X)

  alpha X + Y -> Y
  Add(alpha, X, Y)

  X -> Y
  Copy(X, Y)

  X <-> Y
  Swap(X, Y)

  X.Y
  DotProd(X, Y)
  DotProdConj(X, Y)

  ||X||
  Norm1(X)
  Norm2(X)
  GetMaxAbsIndex(X)

  Omega X
  GenRot(x, y, cos, sin)
  ApplyRot(x, y, cos, sin)

*/


namespace Seldon
{


  /////////
  // MLT //


  template <class T0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const T0 alpha, Vector<T1, Storage1, Allocator1>& X)
  {
    X *= alpha;
  }


  // MLT //
  /////////


  /////////
  // ADD //


  template <class T0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void Add(const T0 alpha,
	   const Vector<T1, Storage1, Allocator1>& X,
	   Vector<T2, Storage2, Allocator2>& Y)
  {
    if (alpha != T0(0))
      {
	T1 alpha_ = alpha;

	int ma = X.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
	CheckDim(X, Y, "Add(alpha, X, Y)");
#endif

	for (int i = 0; i < ma; i++)
	  Y(i) += alpha_ * X(i);
      }
  }


  template <class T0,
            class T1, class Allocator1,
            class T2, class Allocator2>
  void Add(const T0 alpha,
           const Vector<T1, PETScSeq, Allocator1>& X,
           Vector<T2, PETScSeq, Allocator2>& Y)
  {
    if (alpha != T0(0))
      {
        T1 alpha_ = alpha;

        int ma = X.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
        CheckDim(X, Y, "Add(alpha, X, Y)");
#endif

        VecAXPY(Y.GetPetscVector(), alpha_, X.GetPetscVector());
      }
  }


  template <class T0,
            class T1, class Allocator1,
            class T2, class Allocator2>
  void Add(const T0 alpha,
           const Vector<T1, PETScPar, Allocator1>& X,
           Vector<T2, PETScPar, Allocator2>& Y)
  {
    if (alpha != T0(0))
      {
        T1 alpha_ = alpha;

        int ma = X.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
        CheckDim(X, Y, "Add(alpha, X, Y)");
#endif

        VecAXPY(Y.GetPetscVector(), alpha_, X.GetPetscVector());
      }
  }



  template <class T0,
	    class T1, class Allocator1,
	    class T2, class Allocator2>
  void Add(const T0 alpha,
	   const Vector<T1, VectSparse, Allocator1>& X,
	   Vector<T2, VectSparse, Allocator2>& Y)
  {
    if (alpha != T0(0))
      {
	Vector<T1, VectSparse, Allocator1> Xalpha = X;
	Xalpha *= alpha;
	Y.AddInteractionRow(Xalpha.GetSize(),
			    Xalpha.GetIndex(), Xalpha.GetData(), true);
      }
  }


  template <class T0,
	    class T1, class Allocator1,
	    class T2, class Allocator2>
  void Add(const T0 alpha,
	   const Vector<T1, Collection, Allocator1>& X,
	   Vector<T2, Collection, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(X, Y, "Add(X, Y)", "X + Y");
#endif

    for (int i = 0; i < X.GetNvector(); i++)
      Add(alpha, X.GetVector(i), Y.GetVector(i));
  }


  template <class T0,
	    class T1, template <class U1> class Allocator1,
	    class T2, template <class U2> class Allocator2>
  void Add(const T0 alpha,
	   const
           Vector<FloatDouble, DenseSparseCollection, Allocator1<T1> >& X,
	   Vector<FloatDouble, DenseSparseCollection, Allocator2<T2> >& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(X, Y, "Add(X, Y)");
#endif

    Add(alpha, X.GetFloatDense(), Y.GetFloatDense());
    Add(alpha, X.GetFloatSparse(), Y.GetFloatSparse());
    Add(alpha, X.GetDoubleDense(), Y.GetDoubleDense());
    Add(alpha, X.GetDoubleSparse(), Y.GetDoubleSparse());
  }


  template <class T0,
	    class T1, template <class U1> class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void Add(const T0 alpha,
	   const
           Vector<FloatDouble, DenseSparseCollection, Allocator1<T1> >& X,
	   Vector<T2, Storage2, Allocator2>& Y)
  {
    if (alpha != T0(0))
      {
	double alpha_ = alpha;

	int ma = X.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
	CheckDim(X, Y, "Add(alpha, X, Y)");
#endif

	for (int i = 0; i < ma; i++)
	  Y(i) += alpha_ * X(i);

      }
  }


  // ADD //
  /////////


  //////////
  // COPY //


  template <class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void Copy(const Vector<T1, Storage1, Allocator1>& X,
	    Vector<T2, Storage2, Allocator2>& Y)
  {
    Y.Copy(X);
  }


  template <class T1, class Allocator1,
	    class T2, class Allocator2>
  void Copy(const Vector<T1, Collection, Allocator1>& X,
	    Vector<T2, VectFull, Allocator2>& Y)
  {
    Y.Clear();
    for (int i = 0; i < X.GetNvector(); i++)
      Y.PushBack(X.GetVector(i));
  }


  // COPY //
  //////////


  //////////
  // SWAP //


  template <class T, class Storage, class Allocator>
  void Swap(Vector<T, Storage, Allocator>& X,
	    Vector<T, Storage, Allocator>& Y)
  {
    int nx = X.GetM();
    T* data = X.GetData();
    X.Nullify();
    X.SetData(Y.GetM(), Y.GetData());
    Y.Nullify();
    Y.SetData(nx, data);
  }


  template <class T, class Allocator>
  void Swap(Vector<T, VectSparse, Allocator>& X,
	    Vector<T, VectSparse, Allocator>& Y)
  {
    int nx = X.GetM();
    T* data = X.GetData();
    int* index = X.GetIndex();
    X.Nullify();
    X.SetData(Y.GetM(), Y.GetData(), Y.GetIndex());
    Y.Nullify();
    Y.SetData(nx, data, index);
  }


  // SWAP //
  //////////


  /////////////
  // DOTPROD //


  //! Scalar product between two vectors.
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  T1 DotProd(const Vector<T1, Storage1, Allocator1>& X,
	     const Vector<T2, Storage2, Allocator2>& Y)
  {
    T1 value;
    SetComplexZero(value);

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(X, Y, "DotProd(X, Y)");
#endif

    for (int i = 0; i < X.GetM(); i++)
      value += X(i) * Y(i);

    return value;
  }


  //! Scalar product between two vector collections.
  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  typename T1::value_type DotProd(const Vector<T1, Collection, Allocator1>& X,
				  const Vector<T2, Collection, Allocator2>& Y)
  {
    typename T1::value_type value(0);

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(X, Y, "DotProd(X, Y)", "<X, Y>");
#endif

    for (int i = 0; i < X.GetNvector(); i++)
      value += DotProd(X.GetVector(i), Y.GetVector(i));
    return value;
  }


  //! Scalar product between two vector collections.
  template<class T1, template <class U1> class Allocator1,
	   class T2, template <class U2> class Allocator2>
  double
  DotProd(const
          Vector<FloatDouble, DenseSparseCollection, Allocator1<T1> >& X,
          const
          Vector<FloatDouble, DenseSparseCollection, Allocator2<T2> >& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(X, Y, "DotProd(X, Y)");
#endif

    double value(0.);
    value += DotProd(X.GetFloatDense(), Y.GetFloatDense());
    value += DotProd(X.GetFloatSparse(), Y.GetFloatSparse());
    value += DotProd(X.GetDoubleDense(), Y.GetDoubleDense());
    value += DotProd(X.GetDoubleSparse(), Y.GetDoubleSparse());
    return value;
  }


  //! Scalar product between two vectors.
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  T1 DotProdConj(const Vector<T1, Storage1, Allocator1>& X,
		 const Vector<T2, Storage2, Allocator2>& Y)
  {
    return DotProd(X, Y);
  }


  //! Scalar product between two vectors.
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  complex<T1> DotProdConj(const Vector<complex<T1>, Storage1, Allocator1>& X,
			  const Vector<T2, Storage2, Allocator2>& Y)
  {
    complex<T1> value;
    SetComplexZero(value);

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(X, Y, "DotProdConj(X, Y)");
#endif

    for (int i = 0; i < X.GetM(); i++)
      value += conj(X(i)) * Y(i);

    return value;
  }


  //! Scalar product between two sparse vectors.
  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  T1 DotProd(const Vector<T1, VectSparse, Allocator1>& X,
	     const Vector<T2, VectSparse, Allocator2>& Y)
  {
    T1 value;
    SetComplexZero(value);

    int size_x = X.GetSize();
    int size_y = Y.GetSize();
    int kx = 0, ky = 0, pos_x;
    while (kx < size_x)
      {
	pos_x = X.Index(kx);
	while (ky < size_y && Y.Index(ky) < pos_x)
	  ky++;

	if (ky < size_y && Y.Index(ky) == pos_x)
	  value += X.Value(kx) * Y.Value(ky);

	kx++;
      }

    return value;
  }


  //! Scalar product between two sparse vectors.
  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  complex<T1>
  DotProdConj(const Vector<complex<T1>, VectSparse, Allocator1>& X,
	      const Vector<T2, VectSparse, Allocator2>& Y)
  {
    complex<T1> value(0, 0);

    int size_x = X.GetSize();
    int size_y = Y.GetSize();
    int kx = 0, ky = 0, pos_x;
    while (kx < size_x)
      {
	pos_x = X.Index(kx);
	while (ky < size_y && Y.Index(ky) < pos_x)
	  ky++;

	if (ky < size_y && Y.Index(ky) == pos_x)
	  value += conj(X.Value(kx)) * Y.Value(ky);

	kx++;
      }

    return value;
  }


  // DOTPROD //
  /////////////


  ///////////
  // NORM1 //


  template<class T1, class Storage1, class Allocator1>
  T1 Norm1(const Vector<T1, Storage1, Allocator1>& X)
  {
    T1 value(0);

    for (int i = 0; i < X.GetM(); i++)
      value += abs(X(i));

    return value;
  }


  template<class T1, class Storage1, class Allocator1>
  T1 Norm1(const Vector<complex<T1>, Storage1, Allocator1>& X)
  {
    T1 value(0);

    for (int i = 0; i < X.GetM(); i++)
      value += abs(X(i));

    return value;
  }


  template<class T1, class Allocator1>
  T1 Norm1(const Vector<T1, VectSparse, Allocator1>& X)
  {
    T1 value(0);

    for (int i = 0; i < X.GetSize(); i++)
      value += abs(X.Value(i));

    return value;
  }


  template<class T1, class Allocator1>
  T1 Norm1(const Vector<complex<T1>, VectSparse, Allocator1>& X)
  {
    T1 value(0);

    for (int i = 0; i < X.GetSize(); i++)
      value += abs(X.Value(i));

    return value;
  }


  // NORM1 //
  ///////////


  ///////////
  // NORM2 //


  template<class T1, class Storage1, class Allocator1>
  T1 Norm2(const Vector<T1, Storage1, Allocator1>& X)
  {
    T1 value(0);

    for (int i = 0; i < X.GetM(); i++)
      value += X(i) * X(i);

    return sqrt(value);
  }


  template<class T1, class Storage1, class Allocator1>
  T1 Norm2(const Vector<complex<T1>, Storage1, Allocator1>& X)
  {
    T1 value(0);

    for (int i = 0; i < X.GetM(); i++)
      value += real(X(i) * conj(X(i)));

    return sqrt(value);
  }


  template<class T1, class Allocator1>
  T1 Norm2(const Vector<T1, VectSparse, Allocator1>& X)
  {
    T1 value(0);

    for (int i = 0; i < X.GetSize(); i++)
      value += X.Value(i) * X.Value(i);

    return sqrt(value);
  }


  template<class T1, class Allocator1>
  T1 Norm2(const Vector<complex<T1>, VectSparse, Allocator1>& X)
  {
    T1 value(0);

    for (int i = 0; i < X.GetSize(); i++)
      value += real(X.Value(i) * conj(X.Value(i)));

    return sqrt(value);
  }


  // NORM2 //
  ///////////


  ////////////////////
  // GETMAXABSINDEX //


  template<class T, class Storage, class Allocator>
  int GetMaxAbsIndex(const Vector<T, Storage, Allocator>& X)
  {
    return X.GetNormInfIndex();
  }


  // GETMAXABSINDEX //
  ////////////////////


  //////////////
  // APPLYROT //


  //! Computation of rotation between two points.
  template<class T>
  void GenRot(T& a_in, T& b_in, T& c_, T& s_)
  {
    // Old BLAS version.
    T roe;
    if (abs(a_in) > abs(b_in))
      roe = a_in;
    else
      roe = b_in;

    T scal = abs(a_in) + abs(b_in);
    T r, z;
    if (scal != T(0))
      {
	T a_scl = a_in / scal;
	T b_scl = b_in / scal;
	r = scal * sqrt(a_scl * a_scl + b_scl * b_scl);
	if (roe < T(0))
	  r *= T(-1);

	c_ = a_in / r;
	s_ = b_in / r;
	z = T(1);
	if (abs(a_in) > abs(b_in))
	  z = s_;
	else if (abs(b_in) >= abs(a_in) && c_ != T(0))
	  z = T(1) / c_;
      }
    else
      {
	c_ = T(1);
	s_ = T(0);
	r = T(0);
	z = T(0);
      }
    a_in = r;
    b_in = z;
  }


  //! Computation of rotation between two points.
  template<class T>
  void GenRot(complex<T>& a_in, complex<T>& b_in, T& c_, complex<T>& s_)
  {

    T a = abs(a_in), b = abs(b_in);
    if (a == T(0))
      {
	c_ = T(0);
	s_ = complex<T>(1, 0);
	a_in = b_in;
      }
    else
      {
	T scale = a + b;
	T a_scal = abs(a_in / scale);
	T b_scal = abs(b_in / scale);
	T norm = sqrt(a_scal * a_scal + b_scal * b_scal) * scale;

	c_ = a / norm;
	complex<T> alpha = a_in / a;
	s_ = alpha * conj(b_in) / norm;
	a_in = alpha * norm;
      }
    b_in = complex<T>(0, 0);
  }


  //! Rotation of a point in 2-D.
  template<class T>
  void ApplyRot(T& x, T& y, const T c_, const T s_)
  {
    T temp = c_ * x + s_ * y;
    y = c_ * y - s_ * x;
    x = temp;
  }


  //! Rotation of a complex point in 2-D.
  template<class T>
  void ApplyRot(complex<T>& x, complex<T>& y,
		const T& c_, const complex<T>& s_)
  {
    complex<T> temp = s_ * y + c_ * x;
    y = -conj(s_) * x + c_ * y;
    x = temp;
  }


  // APPLYROT //
  //////////////


  //////////////
  // CHECKDIM //


  //! Checks the compatibility of the dimensions.
  /*! Checks that X + Y is possible according to the dimensions of
    the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "X + Y".
  */
  template <class T0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void CheckDim(const Vector<T0, Storage0, Allocator0>& X,
		const Vector<T1, Storage1, Allocator1>& Y,
		string function, string op)
  {
    if (X.GetLength() != Y.GetLength())
      throw WrongDim(function, string("Operation ") + op
		     + string(" not permitted:")
		     + string("\n     X (") + to_str(&X) + string(") is a ")
		     + string("vector of length ") + to_str(X.GetLength())
		     + string(";\n     Y (") + to_str(&Y) + string(") is a ")
		     + string("vector of length ") + to_str(Y.GetLength())
		     + string("."));
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that X + Y is possible according to the dimensions of
    the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "X + Y".
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void CheckDim(const Vector<T0, Vect_Sparse, Allocator0>& X,
		const Vector<T1, Vect_Sparse, Allocator1>& Y,
		string function, string op)
  {
    // The dimension of a Vector<Vect_Sparse> is infinite,
    // so no vector dimension checking has to be done.
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that X + Y is possible according to the dimensions of
    the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "X + Y".
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void CheckDim(const Vector<T0, Collection, Allocator0>& X,
		const Vector<T1, Collection, Allocator1>& Y,
		string function, string op)
  {
    if (X.GetLength() != Y.GetLength())
      throw WrongDim(function, string("Operation ") + op
		     + string(" not permitted:")
		     + string("\n     X (") + to_str(&X) + string(") is a ")
		     + string("vector of length ") + to_str(X.GetLength())
		     + string(";\n     Y (") + to_str(&Y) + string(") is a ")
		     + string("vector of length ") + to_str(Y.GetLength())
		     + string("."));

    if (X.GetNvector() != Y.GetNvector())
      throw WrongDim(function, string("Operation ") + op
		     + string(" not permitted:")
		     + string("\n     X (") + to_str(&X) + string(") is a ")
		     + string("vector of length ") + to_str(X.GetNvector())
		     + string(";\n     Y (") + to_str(&Y) + string(") is a ")
		     + string("vector of length ") + to_str(Y.GetNvector())
		     + string("."));
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that X + Y is possible according to the dimensions of
    the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "X + Y".
  */
  template <class T0, class Allocator0, class Allocator00,
	    class T1, class Allocator1, class Allocator11>
  void CheckDim(const Vector<Vector<T0, Vect_Sparse, Allocator0>,
                Collection, Allocator00>& X,
		const Vector<Vector<T1, Vect_Sparse, Allocator1>,
                Collection, Allocator11>& Y,
		string function, string op)
  {
    if (X.GetNvector() != Y.GetNvector())
      throw WrongDim(function, string("Operation ") + op
		     + string(" not permitted:")
		     + string("\n     X (") + to_str(&X) + string(") is a ")
		     + string("vector of length ") + to_str(X.GetNvector())
		     + string(";\n     Y (") + to_str(&Y) + string(") is a ")
		     + string("vector of length ") + to_str(Y.GetNvector())
		     + string("."));
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that X + Y is possible according to the dimensions of
    the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "X + Y".
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void CheckDim(const Vector<T0, VectFull, Allocator0>& X,
		Vector<T1, Collection, Allocator1>& Y,
		string function, string op)
  {
    if (X.GetLength() != Y.GetM())
      throw WrongDim(function, string("Operation ") + op
		     + string(" not permitted:")
		     + string("\n     X (") + to_str(&X) + string(") is a ")
		     + string("vector of length ") + to_str(X.GetLength())
		     + string(";\n     Y (") + to_str(&Y) + string(") is a ")
		     + string("vector of length ") + to_str(Y.GetM())
		     + string("."));
  }


   //! Checks the compatibility of the dimensions.
  /*! Checks that X + Y is possible according to the dimensions of
    the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "X + Y".
  */
  template <class T0, template <class U0> class Allocator0,
	    class T1, template <class U1> class Allocator1>
  void
  CheckDim(const
           Vector<FloatDouble, DenseSparseCollection, Allocator0<T0> >& X,
           const
           Vector<FloatDouble, DenseSparseCollection, Allocator1<T1> >& Y,
           string function, string op)
  {
    if (X.GetNvector() != Y.GetNvector())
      throw WrongDim(function, string("Operation ") + op
		     + string(" not permitted:")
		     + string("\n     X (") + to_str(&X) + string(") is a ")
		     + string("vector of length ") + to_str(X.GetNvector())
		     + string(";\n     Y (") + to_str(&Y) + string(") is a ")
		     + string("vector of length ") + to_str(Y.GetNvector())
		     + string("."));
  }


  // CHECKDIM //
  //////////////


  ///////////////
  // CONJUGATE //


  //! Sets a vector to its conjugate.
  template<class T, class Prop, class Allocator>
  void Conjugate(Vector<T, Prop, Allocator>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      X(i) = conj(X(i));
  }


  //! Sets a vector to its conjugate.
  template<class T, class Allocator>
  void Conjugate(Vector<T, VectSparse, Allocator>& X)
  {
    for (int i = 0; i < X.GetSize(); i++)
      X.Value(i) = conj(X.Value(i));
  }


  // CONJUGATE //
  ///////////////


} // namespace Seldon.


#endif
