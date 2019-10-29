// Copyright (C) 2012 INRIA
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


#ifndef SELDON_FILE_MATRIX_PETSCMATRIX_HXX


#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"

#include <petscmat.h>


namespace Seldon
{


  //! Matrix class based on PETSc matrix.
  template <class T, class Prop, class Storage,
            class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class PetscMatrix: public Matrix_Base<T, Allocator>
  {
    // Typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef typename Allocator::value_type entry_type;
    typedef typename Allocator::value_type access_type;
    typedef typename Allocator::value_type const_access_type;

    // Attributes.
  protected:
    //! Encapsulated PETSc matrix.
    Mat petsc_matrix_;
    //! The MPI communicator to use.
    MPI_Comm mpi_communicator_;
    //! Boolean to indicate if the inner PETSc matrix is destroyed or not.
    bool petsc_matrix_deallocated_;

    // Methods.
  public:
    // Constructor.
    PetscMatrix();
    PetscMatrix(int i, int j);
    PetscMatrix(Mat& A);

    void SetCommunicator(MPI_Comm mpi_communicator);

    // Destructor.
    ~PetscMatrix();
    void Clear();
    void Nullify();

    Mat& GetPetscMatrix();
    const Mat& GetPetscMatrix() const;

    // Memory management.
    void Resize(int i, int j);
    void SetData(int i, int j, pointer data);

    const_reference Val(int i, int j) const;
    reference Val(int i, int j);
    reference operator[] (int i);
    const_reference operator[] (int i) const;
    void Set(int, int, T);
    void SetBuffer(int, int, T, InsertMode = INSERT_VALUES);
    void Flush() const;
    void GetProcessorRowRange(int& i, int& j) const;
    void Copy(const Mat& A);

    // Convenient functions.
    void Zero();
    void SetIdentity();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
    void FillRand();

    void Print(int a, int b, int m, int n) const;
    void Print(int l) const;

    // Input/output functions.
    void Write(string FileName, bool with_size) const;
    void Write(ostream& FileStream, bool with_size) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;

    void Read(string FileName);
    void Read(istream& FileStream);
    void ReadText(string FileName);
    void ReadText(istream& FileStream);
  };


  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, PETScSeqDense, Allocator>:
    public PetscMatrix<T, Prop, RowMajor, Allocator>
  {
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef typename Allocator::value_type entry_type;
    typedef typename Allocator::value_type access_type;
    typedef typename Allocator::value_type const_access_type;

  public:
    Matrix();
    Matrix(int i, int j);
    Matrix(Mat& A);
    Matrix(const Matrix<T, Prop, PETScSeqDense, Allocator>& A);

    void Reallocate(int i, int j);

    value_type operator() (int i, int j);
    value_type operator() (int i, int j) const;

    void Copy(const Matrix<T, Prop, PETScSeqDense, Allocator>& A);
    Matrix<T, Prop, PETScSeqDense, Allocator>&
    operator= (const Matrix<T, Prop, PETScSeqDense, Allocator>& A);

    void Print() const;
  };


  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, PETScMPIDense, Allocator>:
    public PetscMatrix<T, Prop, RowMajor, Allocator>
  {
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef typename Allocator::value_type entry_type;
    typedef typename Allocator::value_type access_type;
    typedef typename Allocator::value_type const_access_type;

  public:
    Matrix();
    Matrix(int i, int j);
    Matrix(Mat& A);
    Matrix(const Matrix<T, Prop, PETScMPIDense, Allocator>& A);

    void Reallocate(int i, int j, int local_i = PETSC_DECIDE,
                    int local_j = PETSC_DECIDE);

    value_type operator() (int i, int j);
    value_type operator() (int i, int j) const;

    void Copy(const Matrix<T, Prop, PETScMPIDense, Allocator>& A);
    Matrix<T, Prop, PETScMPIDense, Allocator>&
    operator= (const Matrix<T, Prop, PETScMPIDense, Allocator>& A);

    void Print() const;
  };


  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, PETScMPIAIJ, Allocator>:
    public PetscMatrix<T, Prop, RowMajor, Allocator>
  {
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef typename Allocator::value_type entry_type;
    typedef typename Allocator::value_type access_type;
    typedef typename Allocator::value_type const_access_type;

  public:
    Matrix();
    Matrix(int i, int j);
    Matrix(Mat& A);
    Matrix(const Matrix<T, Prop, PETScMPIAIJ, Allocator>& A);

    void Reallocate(int i, int j, int local_i = PETSC_DECIDE,
                    int local_j = PETSC_DECIDE);

    value_type operator() (int i, int j);
    value_type operator() (int i, int j) const;

    void Copy(const Matrix<T, Prop, PETScMPIAIJ, Allocator>& A);
    Matrix<T, Prop, PETScMPIAIJ, Allocator>&
    operator= (const Matrix<T, Prop, PETScMPIAIJ, Allocator>& A);

    template <class T0,  class Allocator0>
    void Copy(const Matrix<T0, General, RowSparse, Allocator0>& A);

    void Print() const;
  };


} // namespace Seldon.


#define SELDON_FILE_MATRIX_PETSCMATRIX_HXX
#endif
