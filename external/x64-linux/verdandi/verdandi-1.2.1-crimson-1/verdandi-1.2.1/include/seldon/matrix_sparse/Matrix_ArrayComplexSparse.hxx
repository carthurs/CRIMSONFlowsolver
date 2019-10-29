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


// To be included by Seldon.hxx

#ifndef SELDON_FILE_MATRIX_ARRAY_COMPLEX_SPARSE_HXX

namespace Seldon
{

  //! Sparse Array-matrix class.
  /*!
    Sparse matrices are defined by: (1) the number of rows and columns;
    (2) the number of non-zero entries; (3) an array of vectors ind
    ind(i) is a vector, which contains indices of columns of the row i;
    (4) an array of vectors val : val(i) is a vector, which contains values of
    the row i
  */
  template <class T, class Prop, class Storage,
	    class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Matrix_ArrayComplexSparse
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef complex<T> entry_type;
    typedef complex<T> access_type;
    typedef complex<T> const_access_type;

    // Attributes.
  protected:
    //! Number of rows.
    int m_;
    //! Number of columns.
    int n_;
    //! real part rows or columns
    Vector<Vector<T, VectSparse, Allocator>, VectFull,
	   NewAlloc<Vector<T, VectSparse, Allocator> > > val_real_;
    //! imaginary part rows or columns
    Vector<Vector<T, VectSparse, Allocator>, VectFull,
	   NewAlloc<Vector<T, VectSparse, Allocator> > > val_imag_;

    // Methods.
  public:
    // Constructors.
    Matrix_ArrayComplexSparse();
    Matrix_ArrayComplexSparse(int i, int j);

    // Destructor.
    ~Matrix_ArrayComplexSparse();
    void Clear();

    // Memory management.
    void Reallocate(int i, int j);
    void Resize(int i, int j);

    // Basic methods.
    int GetM() const;
    int GetN() const;
    int GetM(const SeldonTranspose& status) const;
    int GetN(const SeldonTranspose& status) const;
    int GetRealNonZeros() const;
    int GetImagNonZeros() const;
    int GetRealDataSize() const;
    int GetImagDataSize() const;
    int GetDataSize() const;
    int* GetRealInd(int i) const;
    int* GetImagInd(int i) const;
    T* GetRealData(int i) const;
    T* GetImagData(int i) const;
    Vector<T, VectSparse, Allocator>* GetRealData() const;
    Vector<T, VectSparse, Allocator>* GetImagData() const;

    // Element acess and affectation.
    complex<T> operator() (int i, int j) const;
    complex<T>& Val(int i, int j);
    const complex<T>& Val(int i, int j) const;
    complex<T>& Get(int i, int j);
    const complex<T>& Get(int i, int j) const;

    T& ValReal(int i, int j);
    const T& ValReal(int i, int j) const;
    T& ValImag(int i, int j);
    const T& ValImag(int i, int j) const;
    T& GetReal(int i, int j);
    const T& GetReal(int i, int j) const;
    T& GetImag(int i, int j);
    const T& GetImag(int i, int j) const;

    void Set(int i, int j, const complex<T>& x);

    const T& ValueReal(int num_row,int i) const;
    T& ValueReal(int num_row,int i);
    int IndexReal(int num_row,int i) const;
    int& IndexReal(int num_row,int i);
    const T& ValueImag(int num_row,int i) const;
    T& ValueImag(int num_row,int i);
    int IndexImag(int num_row,int i) const;
    int& IndexImag(int num_row,int i);

    void SetRealData(int, int, Vector<T, VectSparse, Allocator>*);
    void SetImagData(int, int, Vector<T, VectSparse, Allocator>*);
    void SetRealData(int, int, T*, int*);
    void SetImagData(int, int, T*, int*);
    void NullifyReal(int i);
    void NullifyImag(int i);
    void NullifyReal();
    void NullifyImag();

    // Convenient functions.
    void Print() const;
    void Write(string FileName) const;
    void Write(ostream& FileStream) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;
    void Read(string FileName);
    void Read(istream& FileStream);
    void ReadText(string FileName);
    void ReadText(istream& FileStream);

    void Assemble();
    template<class T0>
    void RemoveSmallEntry(const T0& epsilon);

    void SetIdentity();
    void Zero();
    void Fill();
    template <class T0>
    void Fill(const complex<T0>& x);
    template <class T0>
    Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>& operator=
    (const complex<T0>& x);
    void FillRand();

  };


  //! Column-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayColComplexSparse, Allocator> :
    public Matrix_ArrayComplexSparse<T, Prop, ArrayColComplexSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayColComplexSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);

    // Memory management.
    void ClearRealColumn(int i);
    void ClearImagColumn(int i);
    void ReallocateRealColumn(int i, int j);
    void ReallocateImagColumn(int i, int j);
    void ResizeRealColumn(int i, int j);
    void ResizeImagColumn(int i, int j);
    void SwapRealColumn(int i, int i_);
    void SwapImagColumn(int i, int i_);
    void ReplaceRealIndexColumn(int i, IVect& new_index);
    void ReplaceImagIndexColumn(int i, IVect& new_index);

    int GetRealColumnSize(int i) const;
    int GetImagColumnSize(int i) const;
    void PrintRealColumn(int i) const;
    void PrintImagColumn(int i) const;
    void AssembleRealColumn(int i);
    void AssembleImagColumn(int i);

    void AddInteraction(int i, int j, const complex<T>& val);

    template<class Alloc1>
    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<complex<T>, VectFull, Alloc1>& val);
    template<class Alloc1>
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<complex<T>, VectFull,
			      Alloc1>& val);
  };


  //! Row-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowComplexSparse, Allocator> :
    public Matrix_ArrayComplexSparse<T, Prop, ArrayRowComplexSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayRowComplexSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);

    // Memory management.
    void ClearRealRow(int i);
    void ClearImagRow(int i);
    void ReallocateRealRow(int i, int j);
    void ReallocateImagRow(int i, int j);
    void ResizeRealRow(int i, int j);
    void ResizeImagRow(int i, int j);
    void SwapRealRow(int i, int i_);
    void SwapImagRow(int i, int i_);
    void ReplaceRealIndexRow(int i, IVect& new_index);
    void ReplaceImagIndexRow(int i, IVect& new_index);

    int GetRealRowSize(int i) const;
    int GetImagRowSize(int i) const;
    void PrintRealRow(int i) const;
    void PrintImagRow(int i) const;
    void AssembleRealRow(int i);
    void AssembleImagRow(int i);

    void AddInteraction(int i, int j, const complex<T>& val);

    template<class Alloc1>
    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<complex<T>, VectFull, Alloc1>& val);
    template<class Alloc1>
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<complex<T>, VectFull,
			      Alloc1>& val);
  };


  //! Column-major symmetric sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>:
    public Matrix_ArrayComplexSparse<T, Prop, ArrayColSymComplexSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayColSymComplexSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);

    complex<T> operator() (int i, int j) const;

    T& ValReal(int i, int j);
    const T& ValReal(int i, int j) const;
    T& ValImag(int i, int j);
    const T& ValImag(int i, int j) const;
    T& GetReal(int i, int j);
    const T& GetReal(int i, int j) const;
    T& GetImag(int i, int j);
    const T& GetImag(int i, int j) const;

    void Set(int i, int j, const complex<T>& x);

    // Memory management.
    void ClearRealColumn(int i);
    void ClearImagColumn(int i);
    void ReallocateRealColumn(int i, int j);
    void ReallocateImagColumn(int i, int j);
    void ResizeRealColumn(int i, int j);
    void ResizeImagColumn(int i, int j);
    void SwapRealColumn(int i, int i_);
    void SwapImagColumn(int i, int i_);
    void ReplaceRealIndexColumn(int i, IVect& new_index);
    void ReplaceImagIndexColumn(int i, IVect& new_index);

    int GetRealColumnSize(int i) const;
    int GetImagColumnSize(int i) const;
    void PrintRealColumn(int i) const;
    void PrintImagColumn(int i) const;
    void AssembleRealColumn(int i);
    void AssembleImagColumn(int i);

    void AddInteraction(int i, int j, const complex<T>& val);

    template<class Alloc1>
    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<complex<T>, VectFull, Alloc1>& val);
    template<class Alloc1>
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<complex<T>, VectFull,
			      Alloc1>& val);
  };


  //! Row-major symmetric sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>:
    public Matrix_ArrayComplexSparse<T, Prop, ArrayRowSymComplexSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayRowSymComplexSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);

    complex<T> operator() (int i, int j) const;

    T& ValReal(int i, int j);
    const T& ValReal(int i, int j) const;
    T& ValImag(int i, int j);
    const T& ValImag(int i, int j) const;
    T& GetReal(int i, int j);
    const T& GetReal(int i, int j) const;
    T& GetImag(int i, int j);
    const T& GetImag(int i, int j) const;

    void Set(int i, int j, const complex<T>& x);

    // Memory management.
    void ClearRealRow(int i);
    void ClearImagRow(int i);
    void ReallocateRealRow(int i, int j);
    void ReallocateImagRow(int i, int j);
    void ResizeRealRow(int i, int j);
    void ResizeImagRow(int i, int j);
    void SwapRealRow(int i, int i_);
    void SwapImagRow(int i, int i_);
    void ReplaceRealIndexRow(int i, IVect& new_index);
    void ReplaceImagIndexRow(int i, IVect& new_index);

    int GetRealRowSize(int i) const;
    int GetImagRowSize(int i) const;
    void PrintRealRow(int i) const;
    void PrintImagRow(int i) const;
    void AssembleRealRow(int i);
    void AssembleImagRow(int i);

    void AddInteraction(int i, int j, const complex<T>& val);

    template<class Alloc1>
    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<complex<T>, VectFull, Alloc1>& val);
    template<class Alloc1>
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<complex<T>, VectFull,
			      Alloc1>& val);
  };


} // namespace Seldon

#define SELDON_FILE_MATRIX_ARRAY_COMPLEX_SPARSE_HXX
#endif
