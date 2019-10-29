// Copyright (C) 2001-2009 INRIA
// Author(s): Marc Fragu, Vivien Mallet
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


#ifndef SELDON_FILE_VECTOR_HETEROGENEOUSCOLLECTION_HXX


#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"


#ifndef SELDON_DEFAULT_COLLECTION_ALLOCATOR
#define SELDON_DEFAULT_COLLECTION_ALLOCATOR NewAlloc
#endif

namespace Seldon
{

  //! Structure for distributed vectors.
  template <class T, template <class U> class Allocator >
  class Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
    : public Vector_Base<T, Allocator<T> >
  {
    // typedef declarations.
  public:
    typedef Vector<float, VectFull, Allocator<float> > float_dense_v;
    typedef Vector<float, VectSparse, Allocator<float> > float_sparse_v;
    typedef Vector<double, VectFull, Allocator<double> > double_dense_v;
    typedef Vector<double, VectSparse, Allocator<double> > double_sparse_v;

    typedef Vector<float_dense_v, Collection,
		   SELDON_DEFAULT_COLLECTION_ALLOCATOR<float_dense_v> >
    float_dense_c;
    typedef Vector<float_sparse_v, Collection,
		   SELDON_DEFAULT_COLLECTION_ALLOCATOR<float_sparse_v> >
    float_sparse_c;
    typedef Vector<double_dense_v, Collection,
		   SELDON_DEFAULT_COLLECTION_ALLOCATOR<double_dense_v> >
    double_dense_c;
    typedef Vector<double_sparse_v, Collection,
		   SELDON_DEFAULT_COLLECTION_ALLOCATOR<double_sparse_v> >
    double_sparse_c;

    typedef DenseSparseCollection storage;

    // Attributes.
  protected:
    //! Total number of vectors.
    int Nvector_;
    //! For each underlying vectors, index of the corresponding collection.
    Vector<int, VectFull, MallocAlloc<int> > collection_;
    //! Index of the underlying vectors in the inner collection.
    Vector<int, VectFull, MallocAlloc<int> > subvector_;
    //! Lengths of the underlying vectors.
    Vector<int, VectFull, MallocAlloc<int> > length_;
    //! Cumulative sum of the lengths of the underlying vectors.
    Vector<int, VectFull, MallocAlloc<int> > length_sum_;

    //! Pointers of the underlying float dense vectors.
    float_dense_c float_dense_c_;
    //! Pointers of the underlying float sparse vectors.
    float_sparse_c float_sparse_c_;
    //! Pointers of the underlying double dense vectors.
    double_dense_c double_dense_c_;
    //! Pointers of the underlying float sparse vectors.
    double_sparse_c double_sparse_c_;

    //! Indexes of the inner vectors that have a name.
    map<string, int> label_map_;
    //! Names associated with the inner vectors.
    vector<string> label_vector_;

    // Methods.
  public:
    // Constructor.
    explicit Vector();
    Vector(const Vector<FloatDouble, DenseSparseCollection, Allocator<T> >&);

    // Destructor.
    ~Vector();
    void Clear();
    void Deallocate();

    // Management of the vectors.
    void AddVector(const Vector<float, VectFull, Allocator<float> >&);
    void AddVector(const Vector<float, VectSparse, Allocator<float> >&);
    void AddVector(const Vector<double, VectFull, Allocator<double> >&);
    void AddVector(const Vector<double, VectSparse, Allocator<double> >&);

    template <class T0, class Storage0, class Allocator0>
    void AddVector(const Vector<T0, Storage0, Allocator0>&, string name);

    void SetVector(int i,
		   const Vector<float, VectFull, Allocator<float> >&);
    void SetVector(int i,
		   const Vector<float, VectSparse, Allocator<float> >&);
    void SetVector(int i,
		   const Vector<double, VectFull, Allocator<double> >&);
    void SetVector(int i,
		   const Vector<double, VectSparse, Allocator<double> >&);

    template <class T0, class Storage0, class Allocator0>
    void SetVector(int i, const Vector<T0, Storage0, Allocator0>&,
                   string name);
    template <class T0, class Storage0, class Allocator0>
    void SetVector(string name, const Vector<T0, Storage0, Allocator0>&);
    void SetName(int i, string name);

    void Nullify();

    // Basic methods.
    int GetM() const;
    int GetLength() const;
    int GetNvector() const;

    const Vector<int, VectFull, MallocAlloc<int> >& GetVectorLength() const;
    const Vector<int, VectFull, MallocAlloc<int> >& GetLengthSum() const;
    const Vector<int, VectFull, MallocAlloc<int> >& GetCollectionIndex()
      const;
    const Vector<int, VectFull, MallocAlloc<int> >& GetSubvectorIndex() const;

    float_dense_c& GetFloatDense();
    const float_dense_c& GetFloatDense() const;
    float_sparse_c& GetFloatSparse();
    const float_sparse_c& GetFloatSparse() const;
    double_dense_c& GetDoubleDense();
    const double_dense_c& GetDoubleDense() const;
    double_sparse_c& GetDoubleSparse();
    const double_sparse_c& GetDoubleSparse() const;


    void GetVector(int i, float_dense_v& vector) const;
    void GetVector(int i, float_sparse_v& vector) const;
    void GetVector(int i, double_dense_v& vector) const;
    void GetVector(int i, double_sparse_v& vector) const;
    template <class T0, class Storage0, class Allocator0>
    void GetVector(string name, Vector<T0, Storage0, Allocator0>& vector)
      const;

    // Element access and assignment.
    double operator() (int i) const;

    Vector<FloatDouble, DenseSparseCollection, Allocator<T> >& operator=
    (const Vector<FloatDouble, DenseSparseCollection, Allocator<T> >& X);

    void
    Copy(const Vector<FloatDouble, DenseSparseCollection, Allocator<T> >& X);

    template <class T0>
    Vector<FloatDouble, DenseSparseCollection, Allocator<T> >&
    operator*= (const T0& X);

    // Convenient method.
    void Print() const;

    // Input/output functions.
    void Write(string FileName, bool with_size) const;
    void Write(ostream& FileStream, bool with_size) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;

    void Read(string FileName);
    void Read(istream& FileStream);

  protected:
    string GetType(int i) const;
  };

  template <class T, template <class U> class Allocator >
  ostream& operator << (ostream& out,
			const Vector<FloatDouble, DenseSparseCollection,
                        Allocator<T> >& V);

} // namespace Seldon.


#define SELDON_FILE_VECTOR_HETEROGENEOUSCOLLECTION_HXX
#endif
