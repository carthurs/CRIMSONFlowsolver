// Copyright (C) 2011-2012, INRIA
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


#ifndef SELDON_FILE_VECTOR_PETSCVECTOR_HXX

#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"

#include "petscvec.h"


namespace Seldon
{


  //! PETSc vector class.
  /*!
    PETSc vector class.
  */
  template <class T, class Allocator>
  class PETScVector: public Vector_Base<T, Allocator>
  {
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;

    //! Encapsulated PETSc vector.
    Vec petsc_vector_;
    //! The MPI communicator to use.
    MPI_Comm mpi_communicator_;
    //! Boolean to indicate if the inner PETSc vector is destroyed or not.
    bool petsc_vector_deallocated_;

    // Methods.
  public:
    // Constructor.
    explicit PETScVector();
    explicit PETScVector(int i, MPI_Comm mpi_communicator = MPI_COMM_WORLD);
    PETScVector(Vec& petsc_vector);
    PETScVector(const PETScVector<T, Allocator>& A);

    // Destructor.
    ~PETScVector();

    Vec& GetPetscVector();
    const Vec& GetPetscVector() const;
    void SetCommunicator(MPI_Comm mpi_communicator);

    void Clear();

    void Resize(int i);
    void SetData(int i, pointer data);
    void Nullify();

    // Element access.
    value_type operator() (int i) const;
    void SetBuffer(int i, T value, InsertMode insert_mode = INSERT_VALUES);
    void Flush();
    void GetProcessorRange(int& i, int& j) const;
    void Copy(const PETScVector<T, Allocator>& X);
    void Copy(const Vec& petsc_vector);
    void Append(const T& x);

    // Basic functions.
    int GetDataSize();
    int GetLocalM();

    // Convenient functions.
    void Zero();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
    void FillRand();

    // Norms.
    value_type GetNormInf() const;
    int GetNormInfIndex() const;
  };


  template <class T, class Allocator>
  class Vector<T, PETScSeq, Allocator>: public PETScVector<T, Allocator>
  {
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;

    typedef PETScSeq storage;

    // Methods.
  public:
    // Constructor.
    explicit Vector();
    explicit Vector(int i, MPI_Comm mpi_communicator = MPI_COMM_WORLD);
    Vector(Vec& petsc_vector);
    Vector(const Vector<T, PETScSeq, Allocator>& A);
    // Destructor.
    ~Vector();

    void Copy(const Vector<T, PETScSeq, Allocator>& X);
    void Copy(const Vec& petsc_vector);
    // Memory management.
    void Reallocate(int i);

#ifndef SWIG
    Vector<T, PETScSeq, Allocator>& operator= (const Vector<T,
                                               PETScSeq, Allocator>& X);
    template <class T0>
    Vector<T, PETScSeq, Allocator>& operator= (const T0& X);
#endif
    template <class T0>
    Vector<T, PETScSeq, Allocator>& operator*= (const T0& X);

    void Print() const;

    // Input/output functions.
    void Write(string FileName, bool with_size = true) const;
    void Write(ostream& FileStream, bool with_size = true) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;
    void Read(string FileName, bool with_size = true);
    void Read(istream& FileStream, bool with_size = true);
    void ReadText(string FileName);
    void ReadText(istream& FileStream);
  };


#ifndef SWIG
  template <class T, class Allocator>
  ostream& operator << (ostream& out,
                        const Vector<T, PETScSeq, Allocator>& V);
#endif


  template <class T, class Allocator>
  class Vector<T, PETScPar, Allocator>: public PETScVector<T, Allocator>
  {
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;

    typedef PETScPar storage;

    // Methods.
  public:
    // Constructor.
    explicit Vector();
    explicit Vector(int i, MPI_Comm mpi_communicator = MPI_COMM_WORLD);
    explicit Vector(int i, int Nlocal, MPI_Comm mpi_communicator);
    Vector(Vec& petsc_vector);
    Vector(const Vector<T, PETScPar, Allocator>& A);
    // Destructor.
    ~Vector();

    void Copy(const Vector<T, PETScPar, Allocator>& X);
    void Copy(const Vec& petsc_vector);
    // Memory management.
    void Reallocate(int i, int local_size = PETSC_DECIDE);

#ifndef SWIG
    Vector<T, PETScPar, Allocator>& operator= (const Vector<T,
                                               PETScPar, Allocator>& X);
    template <class T0>
    Vector<T, PETScPar, Allocator>& operator= (const T0& X);
#endif
    template <class T0>
    Vector<T, PETScPar, Allocator>& operator*= (const T0& X);

    void Print() const;

    // Input/output functions.
    void Write(string FileName, bool with_size = true) const;
    void Write(ostream& FileStream, bool with_size = true) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;
    void Read(string FileName, bool with_size = true);
    void Read(istream& FileStream, bool with_size = true);
    void ReadText(string FileName);
    void ReadText(istream& FileStream);
  };


#ifndef SWIG
  template <class T, class Allocator>
  ostream& operator << (ostream& out,
                        const Vector<T, PETScPar, Allocator>& V);
#endif



} // namespace Seldon.

#define SELDON_FILE_VECTOR_PETSCVECTOR_HXX
#endif
