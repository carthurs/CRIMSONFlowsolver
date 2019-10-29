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


#ifndef SELDON_FILE_HETEROGENEOUS_MATRIX_COLLECTION_CXX

#include "HeterogeneousMatrixCollection.hxx"

namespace Seldon
{


  ///////////////////////////////////
  // HETEROGENEOUSMATRIXCOLLECTION //
  ///////////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::HeterogeneousMatrixCollection():
    Matrix_Base<double, Allocator<double> >(), Mlocal_(), Mlocal_sum_(1),
    Nlocal_(), Nlocal_sum_(1), collection_(), float_dense_c_(),
    float_sparse_c_(), double_dense_c_(), double_sparse_c_()
  {
    nz_ = 0;
    Mmatrix_ = 0;
    Nmatrix_ = 0;
    Mlocal_sum_.Fill(0);
    Nlocal_sum_.Fill(0);
    collection_.Fill(-1);
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::HeterogeneousMatrixCollection(int i, int j):
    Matrix_Base<double, Allocator<double> >(i, j),
    Mlocal_(i), Mlocal_sum_(i + 1),
    Nlocal_(j), Nlocal_sum_(j + 1), collection_(i, j), float_dense_c_(i, j),
    float_sparse_c_(i, j), double_dense_c_(i, j), double_sparse_c_(i, j)
  {
    nz_ = 0;
    Mmatrix_ = i;
    Nmatrix_ = j;
    Mlocal_.Fill(0);
    Nlocal_.Fill(0);
    Mlocal_sum_.Fill(0);
    Nlocal_sum_.Fill(0);
    collection_.Fill(-1);
  }


  //! Copy constructor.
  /*!
    \param[in] A matrix collection to be copied.
    \note Memory is duplicated: \a A is therefore independent from the current
    instance after the copy.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::HeterogeneousMatrixCollection(const HeterogeneousMatrixCollection<Prop0,
				  Storage0, Prop1, Storage1, Allocator>& A):
    Matrix_Base<double, Allocator<double> >()
  {
    this->Copy(A);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::~HeterogeneousMatrixCollection()
  {
    this->Clear();
  }


  //! Clears the matrix collection without releasing memory.
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Clear()
  {
    float_dense_c_.Nullify();
    float_sparse_c_.Nullify();
    double_dense_c_.Nullify();
    double_sparse_c_.Nullify();

    nz_ = 0;
    Mmatrix_ = 0;
    Nmatrix_ = 0;
    Mlocal_.Clear();
    Nlocal_.Clear();
    Mlocal_sum_.Clear();
    Nlocal_sum_.Clear();
    collection_.Clear();
  }


  //! Clears the matrix collection without releasing memory.
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Nullify()
  {
    float_dense_c_.Nullify();
    float_sparse_c_.Nullify();
    double_dense_c_.Nullify();
    double_sparse_c_.Nullify();

    nz_ = 0;
    Mmatrix_ = 0;
    Nmatrix_ = 0;
    Mlocal_.Clear();
    Nlocal_.Clear();
    Mlocal_sum_.Clear();
    Nlocal_sum_.Clear();
    collection_.Clear();
  }


  //! Clears a given underlying matrix.
  /*!
    \param[in] i row of the underlying matrix to be nullified.
    \param[in] j column of the underlying matrix to be nullified.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Nullify(int i, int j)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("HeterogeneousMatrixCollection::Nullify()",
                     string("Index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("HeterogeneousMatrixCollection::Nullify()",
                     string("Index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif

    switch (collection_(i, j))
      {
      case 0:
	nz_ -= float_dense_c_.GetMatrix(i, j).GetDataSize();
	float_dense_c_.Nullify(i, j);
      case 1:
	nz_ -= float_sparse_c_.GetMatrix(i, j).GetDataSize();
	float_sparse_c_.Nullify(i, j);
	break;
      case 2:
	nz_ -= double_dense_c_.GetMatrix(i, j).GetDataSize();
	double_dense_c_.Nullify(i, j);
	break;
      case 3:
	nz_ -= double_sparse_c_.GetMatrix(i, j).GetDataSize();
	double_sparse_c_.Nullify(i, j);
	break;
      }

    collection_(i, j) = -1;
  }


  //! Deallocates underlying the matrices.
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Deallocate()
  {
    float_dense_c_.Deallocate();
    float_sparse_c_.Deallocate();
    double_dense_c_.Deallocate();
    double_sparse_c_.Deallocate();
    this->~HeterogeneousMatrixCollection();
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of rows.
  /*!
    \return the total number of rows. It is the sum of the number of rows in
    the underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetM() const
  {
    return this->m_;
  }


  //! Returns the number of rows.
  /*!
    \return the total number of rows. It is the sum of the number of rows
    in the underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetMmatrix() const
  {
    return Mmatrix_;
  }


  //! Returns the number of rows in an underlying matrix.
  /*!
    \param[in] i row index of the underlying matrix.
    \return The number of rows in the underlying matrices with row index \a i.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetM(int i) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("HeterogeneousMatrixCollection::GetM()",
                     string("Index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
		     + to_str(i) + ".");
#endif

    return Mlocal_(i);
  }


  //! Returns the number of columns.
  /*!
    \return the total number of columns. It is the sum of the number of
    columns in the underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetN() const
  {
    return this->n_;
  }


  //! Returns the number of columns.
  /*!
    \return the total number of columns. It is the sum of the number of
    columns in the underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetNmatrix() const
  {
    return Nmatrix_;
  }


  //! Returns the number of columns in an underlying matrix.
  /*!
    \param[in] j column index of the underlying matrix.
    \return The number of columns in the underlying matrices with column index
    \a j.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetN(int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("HeterogeneousMatrixCollection::GetN()",
                     string("Index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif

    return Nlocal_(j);
  }


  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetSize() const
  {
    return this->m_ * this->n_;
  }


  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator>
  inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetDataSize() const
  {
    return nz_;
  }


  //! Returns the type of a given underlying matrix.
  /*!
    Type 0 refers to a float dense matrice.
    Type 1 refers to a float sparse matrice.
    Type 2 refers to a double dense matrice.
    Type 3 refers to a double sparse matrice.
    \param[in] i row of the given underlying matrix.
    \param[in] j column of the given underlying matrix.
    \return The type of the underlying matrix.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator>
  inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetType(int i, int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("HeterogeneousMatrixCollection::GetType()",
                     string("Index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("HeterogeneousMatrixCollection::GetType()",
                     string("Index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif
    return collection_(i, j);
  }



  //! Returns the collection of float dense underlying matrices.
  /*!
    \return the collection of float dense underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::float_dense_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetFloatDense()
  {
    return float_dense_c_;
  }


  //! Returns the collection of float dense underlying matrices.
  /*!
    \return the collection of float dense underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline const typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::float_dense_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetFloatDense() const
  {
    return float_dense_c_;
  }


  //! Returns the collection of float sparse underlying matrices.
  /*!
    \return the collection of float sparse underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::float_sparse_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetFloatSparse()
  {
    return float_sparse_c_;
  }


  //! Returns the collection of float sparse underlying matrices.
  /*!
    \return the collection of float sparse underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline const typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::float_sparse_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetFloatSparse() const
  {
    return float_sparse_c_;
  }


  //! Returns the collection of double dense underlying matrices.
  /*!
    \return the collection of double dense underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::double_dense_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetDoubleDense()
  {
    return double_dense_c_;
  }


  //! Returns the collection of double dense underlying matrices.
  /*!
    \return the collection of double dense underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline const typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::double_dense_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetDoubleDense() const
  {
    return double_dense_c_;
  }


  //! Returns the collection of double sparse underlying matrices.
  /*!
    \return the collection of double sparse underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::double_sparse_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetDoubleSparse()
  {
    return double_sparse_c_;
  }


  //! Returns the collection of double sparse underlying matrices.
  /*!
    \return the collection of double sparse underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline const typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::double_sparse_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetDoubleSparse() const
  {
    return double_sparse_c_;
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Reallocates memory to resize the matrix collection.
  /*! On exit, the matrix is a matrix collection with \a i x \a j underlying
    matrices.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
    \warning Depending on your allocator, data may be lost.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Reallocate(int i, int j)
  {
    nz_ = 0;
    Mmatrix_ = i;
    Nmatrix_ = j;
    Mlocal_.Reallocate(i);
    Nlocal_.Reallocate(j);
    Mlocal_sum_.Reallocate(i + 1);
    Nlocal_sum_.Reallocate(j + 1);
    Mlocal_.Fill(0);
    Nlocal_.Fill(0);
    Mlocal_sum_.Fill(0);
    Nlocal_sum_.Fill(0);

    collection_.Reallocate(i, j);
    float_dense_c_.Reallocate(i, j);
    float_sparse_c_.Reallocate(i, j);
    double_dense_c_.Reallocate(i, j);
    double_sparse_c_.Reallocate(i, j);
  }


  //! Sets an underlying  matrix in the matrix collection.
  /*!
    \param[in] i row of the underlying matrix to be set.
    \param[in] j column of the underlying matrix to be set.
    \param[in] matrix new value for the underlying matrix.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::SetMatrix(int i, int j, const typename HeterogeneousMatrixCollection<
	      Prop0, Storage0, Prop1, Storage1, Allocator>::float_dense_m& A)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("HeterogeneousMatrixCollection::"
                     "SetMatrix(float_dense_m)",
                     string("Index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("HeterogeneousMatrixCollection::"
                     "SetMatrix(float_dense_m)",
                     string("Index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
    if ((Mlocal_(i) != 0) && (Mlocal_(i) != A.GetM()))
      throw WrongDim("HeterogeneousMatrixCollection::"
                     "SetMatrix(float_dense_m)",
		     string("The matrix expected should have ")
		     + to_str(this->Mlocal_(i)) + " lines, but has "
		     + to_str(A.GetM()) + " lines.");
    if ((Nlocal_(j) != 0) && (Nlocal_(j) != A.GetN()))
      throw WrongDim("HeterogeneousMatrixCollection::"
                     "SetMatrix(float_dense_m)",
		     string("The matrix expected should have ")
		     + to_str(this->Nlocal_(j)) + " columns, but has "
		     + to_str(A.GetN()) + " columns.");
#endif

    Nullify(i, j);

    collection_(i, j) = 0;

    int Mdiff = A.GetM() - Mlocal_(i);
    int Ndiff = A.GetN() - Nlocal_(j);

    Mlocal_(i) = A.GetM();
    Nlocal_(j) = A.GetN();

    for (int k = i + 1; k < Mmatrix_ + 1; k++)
      Mlocal_sum_(k) += Mdiff;

    for (int k = j + 1; k < Nmatrix_ + 1; k++)
      Nlocal_sum_(k) += Ndiff;

    this->m_ = Mlocal_sum_(Mmatrix_);
    this->n_ = Nlocal_sum_(Nmatrix_);

    float_dense_c_.SetMatrix(i, j, A);
  }


  //! Sets an underlying  matrix in the matrix collection.
  /*!
    \param[in] i row of the underlying matrix to be set.
    \param[in] j column of the underlying matrix to be set.
    \param[in] matrix new value for the underlying matrix.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::SetMatrix(int i, int j, const typename HeterogeneousMatrixCollection<
	      Prop0, Storage0, Prop1, Storage1, Allocator>::float_sparse_m& A)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("HeterogeneousMatrixCollection::"
                     "SetMatrix(float_sparse_m)",
                     string("Index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("HeterogeneousMatrixCollection::"
                     "SetMatrix(float_sparse_m)",
                     string("Index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
    if ((Mlocal_(i) != 0) && (Mlocal_(i) != A.GetM()))
      throw WrongDim("HeterogeneousMatrixCollection::"
                     "SetMatrix(float_sparse_m)",
		     string("The matrix expected should have ")
		     + to_str(this->Mlocal_(i)) + " lines, but has "
		     + to_str(A.GetM()) + " lines.");
    if ((Nlocal_(j) != 0) && (Nlocal_(j) != A.GetN()))
      throw WrongDim("HeterogeneousMatrixCollection::"
                     "SetMatrix(float_sparse_m)",
		     string("The matrix expected should have ")
		     + to_str(this->Nlocal_(j)) + " columns, but has "
		     + to_str(A.GetN()) + " columns.");
#endif

    Nullify(i, j);

    collection_(i, j) = 1;

    int Mdiff = A.GetM() - Mlocal_(i);
    int Ndiff = A.GetN() - Nlocal_(j);

    Mlocal_(i) = A.GetM();
    Nlocal_(j) = A.GetN();

    for (int k = i + 1; k < Mmatrix_ + 1; k++)
      Mlocal_sum_(k) += Mdiff;

    for (int k = j + 1; k < Nmatrix_ + 1; k++)
      Nlocal_sum_(k) += Ndiff;

    this->m_ = Mlocal_sum_(Mmatrix_);
    this->n_ = Nlocal_sum_(Nmatrix_);

    float_sparse_c_.SetMatrix(i, j, A);
  }


  //! Sets an underlying  matrix in the matrix collection.
  /*!
    \param[in] i row of the underlying matrix to be set.
    \param[in] j column of the underlying matrix to be set.
    \param[in] matrix new value for the underlying matrix.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::SetMatrix(int i, int j, const typename HeterogeneousMatrixCollection<
	      Prop0, Storage0, Prop1, Storage1, Allocator>::double_dense_m& A)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("HeterogeneousMatrixCollection::"
                     "SetMatrix(double_dense_m)",
                     string("Index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("HeterogeneousMatrixCollection::"
                     "SetMatrix(double_dense_m)",
                     string("Index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
    if ((Mlocal_(i) != 0) && (Mlocal_(i) != A.GetM()))
      throw WrongDim("HeterogeneousMatrixCollection::"
                     "SetMatrix(double_dense_m)",
		     string("The matrix expected should have ")
		     + to_str(this->Mlocal_(i)) + " lines, but has "
		     + to_str(A.GetM()) + " lines.");
    if ((Nlocal_(j) != 0) && (Nlocal_(j) != A.GetN()))
      throw WrongDim("HeterogeneousMatrixCollection::"
                     "SetMatrix(double_dense_m)",
		     string("The matrix expected should have ")
		     + to_str(this->Nlocal_(j)) + " columns, but has "
		     + to_str(A.GetN()) + " columns.");
#endif

    Nullify(i, j);

    collection_(i, j) = 2;

    int Mdiff = A.GetM() - Mlocal_(i);
    int Ndiff = A.GetN() - Nlocal_(j);

    Mlocal_(i) = A.GetM();
    Nlocal_(j) = A.GetN();

    for (int k = i + 1; k < Mmatrix_ + 1; k++)
      Mlocal_sum_(k) += Mdiff;

    for (int k = j + 1; k < Nmatrix_ + 1; k++)
      Nlocal_sum_(k) += Ndiff;

    this->m_ = Mlocal_sum_(Mmatrix_);
    this->n_ = Nlocal_sum_(Nmatrix_);

    double_dense_c_.SetMatrix(i, j, A);
  }


  //! Sets an underlying  matrix in the matrix collection.
  /*!
    \param[in] i row of the underlying matrix to be set.
    \param[in] j column of the underlying matrix to be set.
    \param[in] matrix new value for the underlying matrix.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::SetMatrix(int i, int j,
              const typename HeterogeneousMatrixCollection< Prop0, Storage0,
              Prop1, Storage1, Allocator>::double_sparse_m& A)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("HeterogeneousMatrixCollection::"
                     "SetMatrix(double_sparse_m)",
                     string("Index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("HeterogeneousMatrixCollection::"
                     "SetMatrix(double_sparse_m)",
                     string("Index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
    if ((Mlocal_(i) != 0) && (Mlocal_(i) != A.GetM()))
      throw WrongDim("HeterogeneousMatrixCollection::"
                     "SetMatrix(double_sparse_m)",
		     string("The matrix expected should have ")
		     + to_str(this->Mlocal_(i)) + " lines, but has "
		     + to_str(A.GetM()) + " lines.");
    if ((Nlocal_(j) != 0) && (Nlocal_(j) != A.GetN()))
      throw WrongDim("HeterogeneousMatrixCollection::"
                     "SetMatrix(double_sparse_m)",
		     string("The matrix expected should have ")
		     + to_str(this->Nlocal_(j)) + " columns, but has "
		     + to_str(A.GetN()) + " columns.");
#endif

    Nullify(i, j);

    collection_(i, j) = 3;

    int Mdiff = A.GetM() - Mlocal_(i);
    int Ndiff = A.GetN() - Nlocal_(j);

    Mlocal_(i) = A.GetM();
    Nlocal_(j) = A.GetN();

    for (int k = i + 1; k < Mmatrix_ + 1; k++)
      Mlocal_sum_(k) += Mdiff;

    for (int k = j + 1; k < Nmatrix_ + 1; k++)
      Nlocal_sum_(k) += Ndiff;

    this->m_ = Mlocal_sum_(Mmatrix_);
    this->n_ = Nlocal_sum_(Nmatrix_);

    double_sparse_c_.SetMatrix(i, j, A);
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access to an underlying matrix.
  /*!
    Returns the underlying matrix (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \param[out] The matrix collection (i, j).
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetMatrix(int i, int j, typename HeterogeneousMatrixCollection<
	      Prop0, Storage0, Prop1, Storage1, Allocator>::float_dense_m& M)
    const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("HeterogeneousMatrixCollection::"
                     "GetMatrix(float_dense_m)",
                     string("Row index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("HeterogeneousMatrixCollection::"
                     "GetMatrix(float_dense_m)",
                     string("Column index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif

    if (collection_(i, j) != 0)
      {
	string matrix_type;
	switch(collection_(i, j))
	  {
	  case 1:
	    matrix_type = "float_sparse";
	    break;
	  case 2:
	    matrix_type = "double_dense";
	    break;
	  case 3:
	    matrix_type = "double_sparse";
	    break;
	  default:
	    throw
              WrongArgument("HeterogeneousMatrixCollection::GetMatrix(i, j,"
                            "Matrix<Float, Dense> M)",
                            "Underlying matrix (" + to_str(i) + " ,"
                            + to_str(j) + " ) not defined.");
	  }

	throw WrongArgument("HeterogeneousMatrixCollection::GetMatrix(i, j, "
			    "Matrix<Float, Dense> M)",
			    string("Wrong type for matrix ")
                            + matrix_type + " M.");
      }

    M.SetData(Mlocal_(i), Nlocal_(j),
              float_dense_c_.GetMatrix(i, j).GetData());
  }


  //! Access to an underlying matrix.
  /*!
    Returns the underlying matrix (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \param[out] The matrix collection (i, j).
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetMatrix(int i, int j, typename HeterogeneousMatrixCollection<
	      Prop0, Storage0, Prop1, Storage1, Allocator>::float_sparse_m& M)
    const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("HeterogeneousMatrixCollection::"
                     "GetMatrix(float_sparse_m)",
                     string("Row index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("HeterogeneousMatrixCollection::"
                     "GetMatrix(float_sparse_m)",
                     string("Column index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif

    if (collection_(i, j) != 1)
      {
	string matrix_type;
	switch(collection_(i, j))
	  {
	  case 0:
	    matrix_type = "float_dense";
	    break;
	  case 2:
	    matrix_type = "double_dense";
	    break;
	  case 3:
	    matrix_type = "double_sparse";
	    break;
	  default:
	    throw
              WrongArgument("HeterogeneousMatrixCollection::GetMatrix(i, j, "
                            "Matrix<Float, Sparse> M)",
                            "Underlying matrix (" + to_str(i) + " ,"
                            + to_str(j) + " ) not defined.");
	  }

	throw WrongArgument("HeterogeneousMatrixCollection::GetMatrix(i, j, "
			    "Matrix<Float, Sparse> M)",
			    string("Wrong type for matrix ")
                            + matrix_type + " M.");
      }

    M.SetData(float_sparse_c_.GetMatrix(i, j).GetM(),
	      float_sparse_c_.GetMatrix(i, j).GetN(),
	      float_sparse_c_.GetMatrix(i, j).GetNonZeros(),
	      float_sparse_c_.GetMatrix(i, j).GetData(),
	      float_sparse_c_.GetMatrix(i, j).GetPtr(),
	      float_sparse_c_.GetMatrix(i, j).GetInd());
  }


  //! Access to an underlying matrix.
  /*!
    Returns the underlying matrix (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \param[out] The matrix collection (i, j).
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetMatrix(int i, int j, typename HeterogeneousMatrixCollection<
	      Prop0, Storage0, Prop1, Storage1, Allocator>::double_dense_m& M)
    const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("HeterogeneousMatrixCollection::"
                     "GetMatrix(double_dense_m)",
                     string("Row index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("HeterogeneousMatrixCollection::"
                     "GetMatrix(double_dense_m)",
                     string("Column index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif

    if (collection_(i, j) != 2)
      {
	string matrix_type;
	switch(collection_(i, j))
	  {
	  case 0:
	    matrix_type = "float_dense";
	    break;
	  case 1:
	    matrix_type = "float_sparse";
	    break;
	  case 3:
	    matrix_type = "double_sparse";
	    break;
	  default:
	    throw
              WrongArgument("HeterogeneousMatrixCollection::GetMatrix(i, j, "
                            "Matrix<Double, Dense> M)",
                            "Underlying matrix (" + to_str(i) + " ,"
                            + to_str(j) + " ) not defined.");
	  }

	throw WrongArgument("HeterogeneousMatrixCollection::GetMatrix(i, j, "
                            "Matrix<Double, Dense> M)",
                            string("Wrong type for matrix ")
                            + matrix_type + " M.");
      }

    M.SetData(Mlocal_(i), Nlocal_(j),
              double_dense_c_.GetMatrix(i, j).GetData());
  }


  //! Access to an underlying matrix.
  /*!
    Returns the underlying matrix (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \param[out] The matrix collection (i, j).
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetMatrix(int i, int j, typename HeterogeneousMatrixCollection<Prop0,
              Storage0, Prop1, Storage1, Allocator>::double_sparse_m& M)
    const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("HeterogeneousMatrixCollection::"
                     "GetMatrix(double_sparse_m)",
                     string("Row index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("HeterogeneousMatrixCollection::"
                     "GetMatrix(double_sparse_m)",
                     string("Column index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif

    if (collection_(i, j) != 3)
      {
	string matrix_type;
	switch(collection_(i, j))
	  {
	  case 0:
	    matrix_type = "float_dense";
	    break;
	  case 1:
	    matrix_type = "float_sparse";
	    break;
	  case 2:
	    matrix_type = "double_dense";
	    break;
	  default:
	    throw
              WrongArgument("HeterogeneousMatrixCollection::GetMatrix(i, j,"
                            "Matrix<Double, Sparse> M)",
                            "Underlying matrix (" + to_str(i) + " ,"
                            + to_str(j) + " ) not defined.");
          }

	throw WrongArgument("HeterogeneousMatrixCollection::GetMatrix(i, j,"
			    "Matrix<Double, Sparse> M)",
			    string("Wrong type for matrix ")
                            + matrix_type + " M.");
      }

    M.Nullify();
    M.SetData(double_sparse_c_.GetMatrix(i, j).GetM(),
	      double_sparse_c_.GetMatrix(i, j).GetN(),
	      double_sparse_c_.GetMatrix(i, j).GetNonZeros(),
	      double_sparse_c_.GetMatrix(i, j).GetData(),
	      double_sparse_c_.GetMatrix(i, j).GetPtr(),
	      double_sparse_c_.GetMatrix(i, j).GetInd());
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \return Element (i, j) of the matrix.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator>
  inline
  double HeterogeneousMatrixCollection<Prop0, Storage0,
                                       Prop1, Storage1, Allocator>
  ::operator() (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->Mlocal_sum_(Mmatrix_))
      throw WrongRow("HeterogeneousMatrixCollection::operator()",
                     string("Index should be in [0, ")
                     + to_str(this->Mlocal_sum_(Mmatrix_) - 1)
                     + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= this->Nlocal_sum_(Nmatrix_))
      throw WrongCol("HeterogeneousMatrixCollection::operator()",
                     string("Index should be in [0, ")
                     + to_str(this->Nlocal_sum_(Nmatrix_) - 1)
                     + "], but is equal to "
                     + to_str(j) + ".");
#endif

    int i_global = 0;
    while (i >= Mlocal_sum_(i_global))
      i_global++;
    i_global--;

    int j_global = 0;
    while (j >= Nlocal_sum_(j_global))
      j_global++;
    j_global--;

    double res = 0.;
    switch(collection_(i_global, j_global))
      {
      case 0:
	res = double(float_dense_c_.GetMatrix(i_global, j_global)
		     (i - Mlocal_sum_(i_global), j - Nlocal_sum_(j_global)));
	break;
      case 1:
	res = double(float_sparse_c_.GetMatrix(i_global, j_global)
		     (i - Mlocal_sum_(i_global), j - Nlocal_sum_(j_global)));
	break;
      case 2:
	res = double_dense_c_.GetMatrix(i_global, j_global)
	  (i - Mlocal_sum_(i_global), j - Nlocal_sum_(j_global));
	break;
      case 3:
	res = double_sparse_c_.GetMatrix(i_global, j_global)
	  (i - Mlocal_sum_(i_global), j - Nlocal_sum_(j_global));
	break;
      default:
	throw
          WrongArgument("HeterogeneousMatrixCollection::operator(int, int)",
                        "Underlying matrix (" + to_str(i) + " ,"
                        + to_str(j) + " ) not defined.");
      }
    return res;
  }


  //! Duplicates a matrix collection (assignment operator).
  /*!
    \param[in] A matrix collection to be copied.
    \note Memory is duplicated: \a A is therefore independent from the current
    instance after the copy.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator>
  inline
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::operator= (const HeterogeneousMatrixCollection<Prop0, Storage0,
               Prop1, Storage1, Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }


  //! Duplicates a matrix collection (assignment operator).
  /*!
    \param[in] A matrix collection to be copied.
    \note Memory is duplicated: \a A is therefore independent from the current
    instance after the copy.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator>
  inline void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Copy(const HeterogeneousMatrixCollection<Prop0, Storage0, Prop1,
         Storage1, Allocator>& A)
  {
    Clear();
    this->nz_ = A.nz_;
    Mmatrix_ = A.Mmatrix_;
    Nmatrix_ = A.Nmatrix_;
    this->m_ = A.GetM();
    this->n_ = A.GetN();

    this->Mlocal_ = A.Mlocal_;
    this->Mlocal_sum_ = A.Mlocal_sum_;
    this->Nlocal_ = A.Nlocal_;
    this->Nlocal_sum_ = A.Nlocal_sum_;

    collection_.Copy(A.collection_);

    float_dense_c_.Reallocate(Mmatrix_, Nmatrix_);
    float_sparse_c_.Reallocate(Mmatrix_, Nmatrix_);
    double_dense_c_.Reallocate(Mmatrix_, Nmatrix_);
    double_sparse_c_.Reallocate(Mmatrix_, Nmatrix_);

    float_dense_m m0a;
    float_sparse_m m1a;
    double_dense_m m2a;
    double_sparse_m m3a;

    for (int i = 0; i < Mmatrix_; i++ )
      for (int j = 0; j < Nmatrix_; j++)
	{
	  switch (A.GetType(i, j))
	    {
	    case 0:
	      A.GetMatrix(i, j, m0a);
	      SetMatrix(i, j, m0a);
	      m0a.Nullify();
	      break;
	    case 1:
	      A.GetMatrix(i, j, m1a);
	      SetMatrix(i, j, m1a);
	      m1a.Nullify();
	      break;
	    case 2:
	      A.GetMatrix(i, j, m2a);
	      SetMatrix(i, j, m2a);
	      m2a.Nullify();
	      break;
	    case 3:
	      A.GetMatrix(i, j, m3a);
	      SetMatrix(i, j, m3a);
	      m3a.Nullify();
	      break;
	    default:
	      throw WrongArgument("Matrix<FloatDouble, DenseSparseCollection>"
				  "::MltAdd(alpha, A, B, beta, C) ",
				  "Underlying matrix  C (" + to_str(i) + " ,"
				  + to_str(j) + " ) not defined.");
	    }
	}
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Displays the matrix collection on the standard output.
  /*!
    Displays elements on the standard output, in text format.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Print() const
  {
    for (int i = 0; i < Mlocal_sum_(Mmatrix_); i++)
      {
	for (int j = 0; j < Nlocal_sum_(Nmatrix_); j++)
	  cout << (*this)(i, j) << endl;
	cout << endl;
      }
  }


  //! Writes the matrix collection in a file.
  /*! Stores the matrix collection in a file in binary format. The number of
    rows of matrices (integer) and the number of columns of matrices (integer)
    are written, and the underlying matrices are then written in the same
    order as in memory (e.g. row-major storage).
    \param[in] FileName output file name.
    \param[in] with_size if set to 'false', the dimensions of the matrix are
    not saved.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Pointers::Write(string FileName)",
                    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream, with_size);

    FileStream.close();
  }


  //! Writes the matrix collection to an output stream.
  /*! Writes the matrix collection to an output stream in binary format.  The
    number of rows of matrices (integer) and the number of columns of matrices
    (integer) are written, and the underlying matrices are then written in the
    same order as in memory (e.g. row-major storage).
    \param[in,out] FileStream output stream.
    \param[in] with_size if set to 'false', the dimensions of the matrix are
    not saved.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Write(ostream& FileStream, bool with_size = true) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("HeterogeneousMatrixCollection"
                    "::Write(ostream& FileStream)",
                    "The stream is not ready.");
#endif

    if (with_size)
      {
        FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&Mmatrix_)),
                         sizeof(int));
        FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&Nmatrix_)),
                         sizeof(int));
      }

    collection_.Write(FileStream, with_size);

    float_dense_m m0a;
    float_sparse_m m1a;
    double_dense_m m2a;
    double_sparse_m m3a;

    int i, j;
    for (i = 0; i < Mmatrix_; i++)
      for (j = 0; j < Nmatrix_; j++)
	{
	  switch (GetType(i, j))
	    {
	    case 0:
	      GetMatrix(i, j, m0a);
	      m0a.Write(FileStream, with_size);
	      m0a.Nullify();
	      break;
	    case 1:
	      throw Undefined("Matrix<FloatDouble, DenseSparseCollection>"
			      "Storage0, Prop1, Storage1, Allocator>"
			      "::Write(ostream& FileStream, bool "
			      "with_size = true) ");
	    case 2:
	      GetMatrix(i, j, m2a);
	      m2a.Write(FileStream, with_size);
	      m2a.Nullify();
	      break;
	    case 3:
	      throw Undefined("Matrix<FloatDouble, DenseSparseCollection>"
			      "Storage0, Prop1, Storage1, Allocator>"
			      "::Write(ostream& FileStream, bool "
			      "with_size = true) ");
	    default:
	      throw WrongArgument("Matrix<FloatDouble, DenseSparseCollection>"
				  "::Write(ostream& FileStream, "
                                  "bool with_size = true) ",
				  "Underlying matrix  A (" + to_str(i) + " ,"
				  + to_str(j) + " ) not defined.");
	    }
	}


#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("HeterogeneousMatrixCollection"
                    "::Write(ostream& FileStream)",
                    "Output operation failed.");
#endif

  }


  //! Writes the matrix collection in a file.
  /*! Stores the matrix in a file in text format. Only the underlying matrices
    are written, without the dimensions. Each row is written on a single line
    and elements of a row are delimited by tabulations.
    \param[in] FileName output file name.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("HeterogeneousMatrixCollection"
                    "::WriteText(string FileName)",
                    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }


  //! Writes the matrix collection to an output stream.
  /*! Stores the matrix to an output stream in text format. Only the
    underlying matrices are written, without the dimensions. Each row is
    written on a single line and elements of a row are delimited by
    tabulations.
    \param[in,out] FileStream output stream.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::WriteText(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("HeterogeneousMatrixCollection"
                    "::WriteText(ostream& FileStream)",
                    "The stream is not ready.");
#endif

    float_dense_m m0a;
    float_sparse_m m1a;
    double_dense_m m2a;
    double_sparse_m m3a;

    int i, j;
    for (i = 0; i < Mmatrix_; i++)
      for (j = 0; j < Nmatrix_; j++)
	{
	  switch (GetType(i, j))
	    {
	    case 0:
	      GetMatrix(i, j, m0a);
	      m0a.WriteText(FileStream);
	      m0a.Nullify();
	      break;
	    case 1:
	      GetMatrix(i, j, m1a);
	      m1a.WriteText(FileStream);
	      m1a.Nullify();
	      break;
	    case 2:
	      GetMatrix(i, j, m2a);
	      m2a.WriteText(FileStream);
	      m2a.Nullify();
	      break;
	    case 3:
	      GetMatrix(i, j, m3a);
	      m3a.WriteText(FileStream);
	      m3a.Nullify();
	      break;
	    default:
	      throw WrongArgument("Matrix<FloatDouble, DenseSparseCollection>"
				  "::Write(ostream& FileStream, "
                                  "bool with_size = true) ",
				  "Underlying matrix  A (" + to_str(i) + " ,"
				  + to_str(j) + " ) not defined.");
	    }
	  FileStream << endl;
	}

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("HeterogeneousMatrixCollection"
                    "::WriteText(ostream& FileStream)",
                    "Output operation failed.");
#endif

  }


  //! Reads the matrix collection from a file.
  /*! Reads a matrix collection stored in binary format in a file.  The number
    of rows of matrices (integer) and the number of columns of matrices
    (integer) are read, and the underlying matrices are then read in the same
    order as it should be in memory (e.g. row-major storage).
    \param[in] FileName input file name.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("HeterogeneousMatrixCollection<Prop0, Storage0, Prop1,"
		    " Storage1, Allocator>::Read(string FileName)",
                    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream);

    FileStream.close();
  }


  //! Reads the matrix collection from an input stream.
  /*! Reads a matrix collection stored in binary format from a stream.  The
    number of rows of matrices (integer) and the number of columns of matrices
    (integer) are read, and the underlying matrices are then read in the same
    order as it should be in memory (e.g. row-major storage).
    \param[in,out] FileStream input stream.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Read(istream& FileStream)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("HeterogeneousMatrixCollection<Prop0, Storage0, Prop1,"
		    " Storage1, Allocator>::Read(istream& FileStream)",
                    "The stream is not ready.");
#endif

    int *new_m, *new_n;
    new_m = new int;
    new_n = new int;

    FileStream.read(reinterpret_cast<char*>(new_m), sizeof(int));
    FileStream.read(reinterpret_cast<char*>(new_n), sizeof(int));

    this->Reallocate(*new_m, *new_n);

    collection_.Read(FileStream);

    float_dense_m m0a;
    float_sparse_m m1a;
    double_dense_m m2a;
    double_sparse_m m3a;
    int i, j;
    for (i = 0; i < Mmatrix_; i++)
      for (j = 0; j < Nmatrix_; j++)
	{
          switch (GetType(i, j))
	    {
	    case 0:
	      m0a.Read(FileStream);
	      SetMatrix(i, j, m0a);
	      m0a.Nullify();
	      break;
	    case 1:
	      throw Undefined("Matrix<FloatDouble, DenseSparseCollection>"
			      "Storage0, Prop1, Storage1, Allocator>"
			      "::Read(istream& FileStream)");
	    case 2:
	      m2a.Read(FileStream);
	      SetMatrix(i, j, m2a);
	      m2a.Nullify();
	      break;
	    case 3:
	      throw Undefined("Matrix<FloatDouble, DenseSparseCollection>"
			      "Storage0, Prop1, Storage1, Allocator>"
			      "::Read(istream& FileStream)");
	      break;
	    default:
	      throw WrongArgument("HeterogeneousMatrixCollection<Prop0, "
				  "Storage0, Prop1, Storage1, Allocator>"
				  "::Read(istream& FileStream) ",
				  "Underlying matrix  A (" + to_str(i) + " ,"
				  + to_str(j) + " ) not defined.");
	    }
	}


    delete new_n;
    delete new_m;

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("HeterogeneousMatrixCollection"
                    "::Read(istream& FileStream)",
                    "Input operation failed.");
#endif

  }


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <template <class U> class Allocator>
  inline
  Matrix<FloatDouble, General, DenseSparseCollection, Allocator<double> >
  ::Matrix():
    HeterogeneousMatrixCollection<General, RowMajor, General,
				  RowSparse, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <template <class U> class Allocator>
  inline
  Matrix<FloatDouble, General, DenseSparseCollection, Allocator<double> >
  ::Matrix(int i, int j):
    HeterogeneousMatrixCollection<General, RowMajor, General,
				  RowSparse, Allocator>(i, j)
  {
  }


} // namespace Seldon.

#define SELDON_FILE_MATRIX_HETEROGENEOUS_COLLECTION_CXX
#endif
