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


#ifndef SELDON_FILE_UMFPACK_CXX

#include "UmfPack.hxx"

namespace Seldon
{

  //! constructor
  template<class T>
  MatrixUmfPack_Base<T>::MatrixUmfPack_Base()
  {
    Symbolic = NULL;
    Numeric = NULL;
    n = 0;

    // allocation of arrays Control and Info
    Control.Reallocate(UMFPACK_CONTROL);
    Info.Reallocate(UMFPACK_INFO);

    print_level = -1;
    transpose = false;
    status_facto = 0;
  }


  //! no message will be displayed by UmfPack
  template<class T>
  void MatrixUmfPack_Base<T>::HideMessages()
  {
    print_level = -1;
    Control(UMFPACK_PRL) = 0;
  }


  //! normal amount of message displayed by UmfPack
  template<class T>
  void MatrixUmfPack_Base<T>::ShowMessages()
  {
    print_level = 1;
    Control(UMFPACK_PRL) = 2;
  }


  template<class T>
  void MatrixUmfPack_Base<T>::ShowFullHistory()
  {
    print_level = 2;
  }


  template<class T>
  int MatrixUmfPack_Base<T>::GetInfoFactorization() const
  {
    return status_facto;
  }


  template<class T>
  void MatrixUmfPack_Base<T>::SelectOrdering(int type)
  {
    Control(UMFPACK_ORDERING) = type;
  }


  template<class T>
  void MatrixUmfPack_Base<T>::SetPermutation(const IVect& permut)
  {
    throw Undefined("MatrixUmfPack_Base::SetPermutation(const Vector&)");
  }


  //! constructor
  MatrixUmfPack<double>::MatrixUmfPack() : MatrixUmfPack_Base<double>()
  {
    ptr_ = NULL;
    ind_ = NULL;
    data_ = NULL;
    umfpack_di_defaults(this->Control.GetData());
    this->HideMessages();
  }


  //! constructor
  MatrixUmfPack<complex<double> >::MatrixUmfPack()
    : MatrixUmfPack_Base<complex<double> >()
  {
    ptr_ = NULL;
    ind_ = NULL;
    data_real_ = NULL;
    data_imag_ = NULL;
    umfpack_zi_defaults(this->Control.GetData());
    this->HideMessages();
  }


  //! destructor
  MatrixUmfPack<double>::~MatrixUmfPack()
  {
    Clear();
  }


  //! we clear present factorization if any
  void MatrixUmfPack<double>::Clear()
  {
    if (this->n > 0)
      {
        // memory used for matrix is released
        free(ptr_);
        free(ind_);
        free(data_);

	// memory for numbering scheme is released
	umfpack_di_free_symbolic(&this->Symbolic) ;

	// memory used by LU factorization is released
	umfpack_di_free_numeric(&this->Numeric) ;

	this->n = 0;
	this->Symbolic = NULL;
	this->Numeric = NULL;
        ptr_ = NULL;
        ind_ = NULL;
        data_ = NULL;
      }
  }


  //! destructor
  MatrixUmfPack<complex<double> >::~MatrixUmfPack()
  {
    Clear();
  }


  //! we clear present factorization if any
  void MatrixUmfPack<complex<double> >::Clear()
  {
    if (this->n > 0)
      {
        // memory used for matrix is released
        free(ptr_);
        free(ind_);
        free(data_real_);
        free(data_imag_);

	// memory for numbering scheme is released
	umfpack_zi_free_symbolic(&this->Symbolic) ;

	// memory used by LU factorization is released
	umfpack_zi_free_numeric(&this->Numeric) ;

        this->n = 0;
	this->Symbolic = NULL;
	this->Numeric = NULL;

        ptr_ = NULL;
        ind_ = NULL;
        data_real_ = NULL;
        data_imag_ = NULL;
      }
  }


  //! factorization of a real matrix in double precision
  template<class Prop, class Storage, class Allocator>
  void MatrixUmfPack<double>::
  FactorizeMatrix(Matrix<double, Prop, Storage, Allocator> & mat,
		  bool keep_matrix)
  {
    // we clear previous factorization
    Clear();

    // conversion to unsymmetric matrix in Column Sparse Column Format
    Matrix<double, General, ColSparse, MallocAlloc<double> > Acsc;
    transpose = false;

    this->n = mat.GetM();
    Copy(mat, Acsc);
    if (!keep_matrix)
      mat.Clear();

    // we retrieve pointers of Acsc and nullify this object
    ptr_ = Acsc.GetPtr();
    ind_ = Acsc.GetInd();
    data_ = Acsc.GetData();
    Acsc.Nullify();

    // factorization with UmfPack
    umfpack_di_symbolic(this->n, this->n, ptr_, ind_, data_, &this->Symbolic,
                        this->Control.GetData(), this->Info.GetData());

    // we display informations about the performed operation
    status_facto =
      umfpack_di_numeric(ptr_, ind_, data_,
                         this->Symbolic, &this->Numeric,
                         this->Control.GetData(), this->Info.GetData());

    // we display informations about the performed operation
    if (print_level > 1)
      {
	umfpack_di_report_status(this->Control.GetData(), status_facto);
	umfpack_di_report_info(this->Control.GetData(),this->Info.GetData());
      }

    if (print_level > 0)
      {
	int size_mem = (this->Info(UMFPACK_SYMBOLIC_SIZE)
			+ this->Info(UMFPACK_NUMERIC_SIZE_ESTIMATE))
	  *this->Info(UMFPACK_SIZE_OF_UNIT);

	cout << "Memory used to store LU factors: "
	     << double(size_mem)/(1024*1024) << " MB" << endl;
      }
  }


  //! Symbolic factorization
  template<class Prop, class Allocator>
  void MatrixUmfPack<double>
  ::PerformAnalysis(Matrix<double, Prop, RowSparse, Allocator> & mat)
  {
    // we clear previous factorization
    Clear();

    // conversion to unsymmetric matrix in Column Sparse Row Format
    Matrix<double, General, RowSparse, MallocAlloc<double> > Acsr;
    transpose = true;

    Copy(mat, Acsr);

    // we retrieve pointers of Acsc and nullify this object
    this->n = mat.GetM();
    ptr_ = Acsr.GetPtr();
    ind_ = Acsr.GetInd();
    data_ = Acsr.GetData();
    Acsr.Nullify();

    // factorization with UmfPack
    umfpack_di_symbolic(this->n, this->n, ptr_, ind_, data_, &this->Symbolic,
			this->Control.GetData(), this->Info.GetData());

  }


  //! Numerical factorization
  template<class Prop, class Allocator>
  void MatrixUmfPack<double>
  ::PerformFactorization(Matrix<double, Prop, RowSparse, Allocator> & mat)
  {
    // we copy values
    double* data = mat.GetData();
    for (int i = 0; i < mat.GetDataSize(); i++)
      data_[i] = data[i];

    status_facto =
      umfpack_di_numeric(ptr_, ind_, data_,
			 this->Symbolic, &this->Numeric,
			 this->Control.GetData(), this->Info.GetData());

    // we display informations about the performed operation
    if (print_level > 1)
      {
	umfpack_di_report_status(this->Control.GetData(), status_facto);
	umfpack_di_report_info(this->Control.GetData(),this->Info.GetData());
      }
  }


  //! resolution of A y = x (result is overwritten in x)
  template<class Allocator2>
  void MatrixUmfPack<double>::Solve(Vector<double, VectFull, Allocator2>& x)
  {
    // local copy of x
    Vector<double, VectFull, Allocator2> b(x);

    int sys = UMFPACK_A;
    if (transpose)
      sys = UMFPACK_Aat;

    int status
      = umfpack_di_solve(sys, ptr_, ind_, data_, x.GetData(),
			 b.GetData(), this->Numeric, this->Control.GetData(),
			 this->Info.GetData());

    // we display informations about the performed operation
    if (print_level > 1)
      umfpack_di_report_status(this->Control.GetData(), status);
  }


  template<class StatusTrans, class Allocator2>
  void MatrixUmfPack<double>::Solve(const StatusTrans& TransA,
				    Vector<double, VectFull, Allocator2>& x)
  {
    if (TransA.NoTrans())
      {
	Solve(x);
	return;
      }

    // local copy of x
    Vector<double, VectFull, Allocator2> b(x);

    int sys = UMFPACK_Aat;
    if (transpose)
      sys = UMFPACK_A;

    int status
      = umfpack_di_solve(sys, ptr_, ind_, data_, x.GetData(),
			 b.GetData(), this->Numeric, this->Control.GetData(),
			 this->Info.GetData());

    // We display information about the performed operation.
    if (print_level > 1)
      umfpack_di_report_status(this->Control.GetData(), status);
  }


  //! LU factorization using UmfPack in double complex precision
  template<class Prop, class Storage,class Allocator>
  void MatrixUmfPack<complex<double> >::
  FactorizeMatrix(Matrix<complex<double>, Prop, Storage, Allocator> & mat,
		  bool keep_matrix)
  {
    Clear();

    this->n = mat.GetM();
    // conversion to CSC format
    Matrix<complex<double>, General, ColSparse,
      MallocAlloc<complex<double> > > Acsc;
    transpose = false;

    Copy(mat, Acsc);
    if (!keep_matrix)
      mat.Clear();

    int nnz = Acsc.GetDataSize();
    complex<double>* data = Acsc.GetData();
    int* ptr = Acsc.GetPtr();
    int* ind = Acsc.GetInd();
    Vector<double, VectFull, MallocAlloc<double> >
      ValuesReal(nnz), ValuesImag(nnz);

    Vector<int, VectFull, MallocAlloc<int> > Ptr(this->n+1), Ind(nnz);

    for (int i = 0; i < nnz; i++)
      {
	ValuesReal(i) = real(data[i]);
	ValuesImag(i) = imag(data[i]);
	Ind(i) = ind[i];
      }

    for (int i = 0; i <= this->n; i++)
      Ptr(i) = ptr[i];

    // we clear intermediary matrix Acsc
    Acsc.Clear();

    // retrieve pointers and nullify Seldon vectors
    data_real_ = ValuesReal.GetData();
    data_imag_ = ValuesImag.GetData();
    ptr_ = Ptr.GetData();
    ind_ = Ind.GetData();
    ValuesReal.Nullify(); ValuesImag.Nullify();
    Ptr.Nullify(); Ind.Nullify();

    // we call UmfPack
    umfpack_zi_symbolic(this->n, this->n, ptr_, ind_,
			data_real_, data_imag_,
			&this->Symbolic, this->Control.GetData(),
			this->Info.GetData());

    status_facto
      = umfpack_zi_numeric(ptr_, ind_, data_real_, data_imag_,
			   this->Symbolic, &this->Numeric,
			   this->Control.GetData(), this->Info.GetData());

    if (print_level > 1)
      {
	umfpack_zi_report_status(this->Control.GetData(), status_facto);
	umfpack_zi_report_info(this->Control.GetData(), this->Info.GetData());
      }

    if (print_level > 0)
      {
	int size_mem = (this->Info(UMFPACK_SYMBOLIC_SIZE)
			+ this->Info(UMFPACK_NUMERIC_SIZE_ESTIMATE))
	  *this->Info(UMFPACK_SIZE_OF_UNIT);

	cout << "Estimated memory used to store LU factors: "
	     << double(size_mem)/(1024*1024) << " MB" << endl;
      }
  }


  //! solves linear system in complex double precision using UmfPack
  template<class Allocator2>
  void MatrixUmfPack<complex<double> >::
  Solve(Vector<complex<double>, VectFull, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  //! Solves linear system in complex double precision using UmfPack.
  template<class StatusTrans, class Allocator2>
  void MatrixUmfPack<complex<double> >::
  Solve(const StatusTrans& TransA,
	Vector<complex<double>, VectFull, Allocator2>& x)
  {
    int m = x.GetM();
    // creation of vectors
    Vector<double> b_real(m), b_imag(m);
    for (int i = 0; i < m; i++)
      {
	b_real(i) = real(x(i));
	b_imag(i) = imag(x(i));
      }

    Vector<double> x_real(m), x_imag(m);
    x_real.Zero();
    x_imag.Zero();

    int sys = UMFPACK_A;
    if (TransA.Trans())
      sys = UMFPACK_Aat;

    int status
      = umfpack_zi_solve(sys, ptr_, ind_, data_real_, data_imag_,
			 x_real.GetData(), x_imag.GetData(),
			 b_real.GetData(), b_imag.GetData(),
			 this->Numeric,
			 this->Control.GetData(), this->Info.GetData());
    b_real.Clear();
    b_imag.Clear();

    for (int i = 0; i < m; i++)
      x(i) = complex<double>(x_real(i), x_imag(i));

    if (print_level > 1)
      umfpack_zi_report_status(this->Control.GetData(), status);
  }


  //! LU factorization
  template<class T, class Prop, class Storage, class Allocator>
  void GetLU(Matrix<T, Prop, Storage, Allocator>& A, MatrixUmfPack<T>& mat_lu,
	     bool keep_matrix = false)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }


  //! LU resolution
  template<class T, class Allocator>
  void SolveLU(MatrixUmfPack<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }


  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
               MatrixUmfPack<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }

}

#define SELDON_FILE_UMFPACK_CXX
#endif
