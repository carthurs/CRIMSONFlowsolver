// Copyright (C) 2003-2009 Marc Duruflé
// Copyright (C) 2001-2009 Vivien Mallet
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


// To be included by Seldon.hxx

#ifndef SELDON_FILE_MATRIX_ARRAY_SPARSE_HXX

namespace Seldon
{

  //! Sparse Array-matrix class.
  /*!
    Sparse matrices are defined by: (1) the number of rows and columns;
    (2) the number of non-zero entries; (3) an array of vectors ind
    ind(i) is a vector, which contains indices of columns of the row i
    (4) an array of vectors val : val(i) is a vector, which contains values of
    the row i
  */
  template <class T, class Prop, class Storage,
	    class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Matrix_ArraySparse
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef T entry_type;
    typedef T& access_type;
    typedef T const_access_type;

    // Attributes.
  protected:
    //! Number of rows.
    int m_;
    //! Number of columns.
    int n_;
    //! rows or columns
    Vector<Vector<T, VectSparse, Allocator>, VectFull,
	   NewAlloc<Vector<T, VectSparse, Allocator> > > val_;

  public:
    // Constructors.
    Matrix_ArraySparse();
    Matrix_ArraySparse(int i, int j);

    // Destructor.
    ~Matrix_ArraySparse();
    void Clear();

    // Memory management.
    void Reallocate(int i, int j);
    void Resize(int i, int j);

    // Basic methods.
    int GetM() const;
    int GetN() const;
    int GetM(const SeldonTranspose& status) const;
    int GetN(const SeldonTranspose& status) const;
    int GetNonZeros() const;
    int GetDataSize() const;
    int* GetIndex(int i) const;
    T* GetData(int i) const;

    Vector<T, VectSparse, Allocator>* GetData() const;

    // Element acess and affectation.
    T operator() (int i, int j) const;
    T& Get(int i, int j);
    const T& Get(int i, int j) const;
    T& Val(int i, int j);
    const T& Val(int i, int j) const;
    void Set(int i, int j, const T& x);

    const T& Value(int num_row, int i) const;
    T& Value(int num_row, int i);
    int Index(int num_row, int i) const;
    int& Index(int num_row, int i);

    void SetData(int, int, Vector<T, VectSparse, Allocator>*);
    void SetData(int, int, T*, int*);
    void Nullify(int i);
    void Nullify();

    // Convenient functions.
    void Print() const;
    void Assemble();
    template<class T0>
    void RemoveSmallEntry(const T0& epsilon);

    void SetIdentity();
    void Zero();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
    template <class T0>
    Matrix_ArraySparse<T, Prop, Storage, Allocator>& operator= (const T0& x);
    void FillRand();

    // Input/output functions.
    void Write(string FileName) const;
    void Write(ostream& FileStream) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;
    void Read(string FileName);
    void Read(istream& FileStream);
    void ReadText(string FileName);
    void ReadText(istream& FileStream);

  };


  //! Column-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayColSparse, Allocator> :
    public Matrix_ArraySparse<T, Prop, ArrayColSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayColSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);

    // Memory management.
    void ClearColumn(int i);
    void ReallocateColumn(int i, int j);
    void ResizeColumn(int i, int j);
    void SwapColumn(int i, int i_);
    void ReplaceIndexColumn(int i, IVect& new_index);

    int GetColumnSize(int i) const;
    void PrintColumn(int i) const;
    void AssembleColumn(int i);

    void AddInteraction(int i, int j, const T& val);

    void AddInteractionRow(int, int, int*, T*);
    void AddInteractionColumn(int, int, int*, T*);

    template<class Alloc1>
    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<T, VectFull, Alloc1>& val);
    template<class Alloc1>
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<T, VectFull, Alloc1>& val);
  };


  //! Row-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowSparse, Allocator> :
    public Matrix_ArraySparse<T, Prop, ArrayRowSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayRowSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);

    // Memory management.
    void ClearRow(int i);
    void ReallocateRow(int i, int j);
    void ResizeRow(int i, int j);
    void SwapRow(int i, int i_);
    void ReplaceIndexRow(int i, IVect& new_index);

    int GetRowSize(int i) const;
    void PrintRow(int i) const;
    void AssembleRow(int i);

    void AddInteraction(int i, int j, const T& val);

    void AddInteractionRow(int, int, int*, T*);
    void AddInteractionColumn(int, int, int*, T*);

    template<class Alloc1>
    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<T, VectFull, Alloc1>& val);
    template<class Alloc1>
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<T, VectFull, Alloc1>& val);
  };

  //! Column-major symmetric sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayColSymSparse, Allocator>:
    public Matrix_ArraySparse<T, Prop, ArrayColSymSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayColSymSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);

    T operator() (int i, int j) const;
    T& Get(int i, int j);
    const T& Get(int i, int j) const;
    T& Val(int i, int j);
    const T& Val(int i, int j) const;
    void Set(int i, int j, const T& x);

    // Memory management.
    void ClearColumn(int i);
    void ReallocateColumn(int i, int j);
    void ResizeColumn(int i, int j);
    void SwapColumn(int i, int i_);
    void ReplaceIndexColumn(int i, IVect& new_index);

    int GetColumnSize(int i) const;
    void PrintColumn(int i) const;
    void AssembleColumn(int i);

    void AddInteraction(int i, int j, const T& val);

    void AddInteractionRow(int, int, int*, T*);
    void AddInteractionColumn(int, int, int*, T*);

    template<class Alloc1>
    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<T, VectFull, Alloc1>& val);
    template<class Alloc1>
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<T, VectFull, Alloc1>& val);
  };


  //! Row-major symmetric sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowSymSparse, Allocator>:
    public Matrix_ArraySparse<T, Prop, ArrayRowSymSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayRowSymSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);

    T operator() (int i, int j) const;
    T& Get(int i, int j);
    const T& Get(int i, int j) const;
    T& Val(int i, int j);
    const T& Val(int i, int j) const;
    void Set(int i, int j, const T& x);

    // Memory management.
    void ClearRow(int i);
    void ReallocateRow(int i, int j);
    void ResizeRow(int i, int j);
    void SwapRow(int i, int i_);
    void ReplaceIndexRow(int i, IVect& new_index);

    int GetRowSize(int i) const;
    void PrintRow(int i) const;
    void AssembleRow(int i);

    void AddInteraction(int i, int j, const T& val);

    void AddInteractionRow(int, int, int*, T*);
    void AddInteractionColumn(int, int, int*, T*);

    template<class Alloc1>
    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<T, VectFull, Alloc1>& val);
    template<class Alloc1>
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<T, VectFull, Alloc1>& val);
  };

} // namespace Seldon

#define SELDON_FILE_MATRIX_ARRAY_SPARSE_HXX
#endif
