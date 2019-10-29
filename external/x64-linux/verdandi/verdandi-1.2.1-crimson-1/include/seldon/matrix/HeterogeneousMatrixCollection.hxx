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


#ifndef SELDON_FILE_HETEROGENEOUS_MATRIX_COLLECTION_HXX


#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"

namespace Seldon
{


  //! Matrix class made of an heterogeneous collection of matrices.
  /*!
    A collection can refer to matrices of different types : float, double,
    dense and sparse.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator>
  class HeterogeneousMatrixCollection:
    public Matrix_Base<double, Allocator<double> >
  {

    // Typedef declarations.
  public:
    typedef Matrix<float, Prop0, Storage0, Allocator<float> > float_dense_m;
    typedef Matrix<float, Prop1, Storage1, Allocator<float> > float_sparse_m;
    typedef Matrix<double, Prop0, Storage0, Allocator<double> >
    double_dense_m;
    typedef Matrix<double, Prop1, Storage1, Allocator<double> >
    double_sparse_m;

    typedef Matrix<float_dense_m, General, RowMajorCollection,
		   NewAlloc<float_dense_m> > float_dense_c;
    typedef Matrix<float_sparse_m, General, RowMajorCollection,
		   NewAlloc<float_sparse_m> > float_sparse_c;
    typedef Matrix<double_dense_m, General, RowMajorCollection,
		   NewAlloc<double_dense_m> > double_dense_c;
    typedef Matrix<double_sparse_m, General, RowMajorCollection,
		   NewAlloc<double_sparse_m> > double_sparse_c;

    // Attributes.
  protected:
    //! Number of non-zero elements.
    int nz_;
    //! Number of rows of matrices.
    int Mmatrix_;
    //! Number of columns of matrices.
    int Nmatrix_;
    //! Number of rows in the underlying matrices.
    Vector<int, VectFull, CallocAlloc<int> > Mlocal_;
    //! Cumulative number of rows in the underlying matrices.
    Vector<int, VectFull, CallocAlloc<int> > Mlocal_sum_;
    //! Number of columns in the underlying matrices.
    Vector<int, VectFull, CallocAlloc<int> > Nlocal_;
    //! Cumulative number of columns in the underlying matrices.
    Vector<int, VectFull, CallocAlloc<int> > Nlocal_sum_;

    //! Type of the underlying matrices.
    /*!
      Type 0 refers to float dense matrices.
      Type 1 refers to float sparse matrices.
      Type 2 refers to double dense matrices.
      Type 3 refers to double sparse matrices.
    */
    Matrix<int, General, RowMajor, CallocAlloc<int> > collection_;

    //! Pointers of the underlying float dense matrices.
    float_dense_c float_dense_c_;
    //! Pointers of the underlying float sparse matrices.
    float_sparse_c float_sparse_c_;
    //! Pointers of the underlying double dense matrices.
    double_dense_c double_dense_c_;
    //! Pointers of the underlying double sparse matrices.
    double_sparse_c double_sparse_c_;


    // Methods.
  public:
    // Constructor.
    HeterogeneousMatrixCollection();
    HeterogeneousMatrixCollection(int i, int j);
    HeterogeneousMatrixCollection
    (const HeterogeneousMatrixCollection<Prop0, Storage0, Prop1,
     Storage1, Allocator>& A);

    // Destructor.
    ~HeterogeneousMatrixCollection();
    void Clear();
    void Nullify();
    void Nullify(int i, int j);
    void Deallocate();

    // Basic methods.
    int GetM() const;
    int GetMmatrix() const;
    int GetM(int i) const;
    int GetN() const;
    int GetNmatrix() const;
    int GetN(int j) const;
    int GetSize() const;
    int GetDataSize() const;
    int GetType(int i, int j) const;

    float_dense_c& GetFloatDense();
    const float_dense_c& GetFloatDense() const;
    float_sparse_c& GetFloatSparse();
    const float_sparse_c& GetFloatSparse() const;
    double_dense_c& GetDoubleDense();
    const double_dense_c& GetDoubleDense() const;
    double_sparse_c& GetDoubleSparse();
    const double_sparse_c& GetDoubleSparse() const;

    // Memory management.
    void Reallocate(int i, int j);

    // Management of the matrices.
    void SetMatrix(int m, int n, const float_dense_m&);
    void SetMatrix(int m, int n, const float_sparse_m&);
    void SetMatrix(int m, int n, const double_dense_m&);
    void SetMatrix(int m, int n, const double_sparse_m&);

    // Element access and affectation.
    void GetMatrix(int m, int n, float_dense_m&) const;
    void GetMatrix(int m, int n, float_sparse_m&) const;
    void GetMatrix(int m, int n, double_dense_m&) const;
    void GetMatrix(int m, int n, double_sparse_m&) const;


    double operator() (int i, int j) const;

    HeterogeneousMatrixCollection<Prop0, Storage0, Prop1,
				  Storage1, Allocator>&
    operator= (const HeterogeneousMatrixCollection<Prop0, Storage0,
	       Prop1, Storage1, Allocator>& A);

    void Copy(const HeterogeneousMatrixCollection<Prop0, Storage0,
	      Prop1, Storage1, Allocator>& A);

    // Convenient functions.
    void Print() const;

    // Input/output functions.
    void Write(string FileName, bool with_size) const;
    void Write(ostream& FileStream, bool with_size) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;

    void Read(string FileName);
    void Read(istream& FileStream);

  };


  //! Heterogeneous matrix collection class.
  template <template <class U> class Allocator>
  class Matrix<FloatDouble, General,
               DenseSparseCollection, Allocator<double> >:
    public HeterogeneousMatrixCollection<General, RowMajor, General,
					 RowSparse, Allocator >
  {
    // typedef declaration.
  public:
    typedef General property;
    typedef DenseSparseCollection storage;
    typedef Allocator<double> allocator;

  public:
    Matrix();
    Matrix(int i, int j);
  };


} // namespace Seldon.


#define SELDON_FILE_HETEROGENEOUS_MATRIX_COLLECTION_HXX
#endif
