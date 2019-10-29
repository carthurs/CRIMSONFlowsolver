// Copyright (C) 2008, INRIA
// Author(s): Vivien Mallet
//
// This file is part of the data assimilation library Verdandi.
//
// Verdandi is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Verdandi. If not, see http://www.gnu.org/licenses/.
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_FILE_SHARE_USEFULFUNCTION_HXX

#include <sstream>
#include <string>
#include <iostream>

namespace Verdandi
{


    template <class T, class TM>
    T interpolate(T x_min, T Delta_x, T y_min, T Delta_y,
                  const Matrix<TM>& input, T x, T y);


    void get_position(int index, const Vector<int>& shape,
                      Vector<int>& position);
    template <class T>
    void get_coordinate(int index, const Vector<T>& minimum,
                        const Vector<T>& step, const Vector<int>& shape,
                        Vector<T>& coordinate);

    /*** From Talos library ***/

    bool is_num(const string& s);
    bool is_integer(const string& s);
    bool is_unsigned_integer(const string& s);
    string trim(string str, string delimiters = " \n\t");
    template <class T>
    void split(string str, vector<T>& vect, string delimiters = " \n\t");
    vector<string> split(string str, string delimiters = " \n\t");
    string find_replace(string str, string old_str, string new_str);
    template<class T>
    bool is_equal(T x, T y, T epsilon = 1.e-6);
    template<class T>
    bool is_multiple(T x, T d, T epsilon = 1.e-6);
    string upper_case(string str);

    /*** Builds a diagonal sparse matrix ***/

    template <class T>
    void build_diagonal_sparse_matrix(int size, T diagonal_value,
                                      Matrix<T, General, RowSparse>& matrix);

    /*** Linear algebra ***/

    template <class T, class Allocator>
    void GetCholesky(Matrix<T, General, RowMajor, Allocator>& A);

    template <class T, class Allocator>
    void GetCholesky(Matrix<T, General, RowSparse, Allocator>& A);

    template <class T, class Allocator>
    void GetInverse(Matrix<T, General, RowSparse, Allocator>& A);

    template <class T, class Prop0, class Allocator0, class Allocator1>
    void GetAndSolveLU(Matrix<T, Prop0, ColSparse, Allocator0>& M,
                       Vector<T, VectFull, Allocator1>& Y);

    template <class T, class Prop0, class Allocator0, class Allocator1>
    void GetAndSolveLU(Matrix<T, Prop0, RowSparse, Allocator0>& M,
                       Vector<T, VectFull, Allocator1>& Y);

    template <class T0, class Allocator0,
              class T1, class Allocator1>
    void Copy(const Matrix<T0, General, RowMajor, Allocator0>& A,
              Matrix<T1, General, RowSymPacked, Allocator1>& B);
    template <class T, class Allocator>
    void Copy(const Matrix<T, General, RowSparse, Allocator>& A,
              Matrix<T, General, RowMajor, Allocator>& A_dense);

    template <class T0,
              class T1, class Allocator1,
              class T2, class Allocator2>
    void Add(const T0 alpha, const Vector<T1, Collection, Allocator1>& X,
             Vector<T2, VectFull, Allocator2>& Y);

    template <class T0, class Allocator0, class T1, class Allocator1>
    void SetCol(const Vector<T1, VectFull, Allocator1>& X,
                int i, Matrix<T0, General, RowSparse, Allocator0>& M);

    template <class T, class Allocator>
    void ConvertDenseToArrayRowSparse(
        const Matrix<T, General, RowMajor, Allocator>& A_dense,
        Matrix<T, General, ArrayRowSparse, Allocator>& A_array);
    template <class T, class Allocator>
    void ConvertArrayRowSparseToDense(
        const Matrix<T, General, ArrayRowSparse, Allocator>& A_array,
        Matrix<T, General, RowMajor, Allocator>& A_dense);
    template <class T, class Allocator>
    void ConvertRowSparseToArrayRowSparse(
        const Matrix<T, General, RowSparse, Allocator>& A,
        Matrix<T, General, ArrayRowSparse, Allocator>& A_dense);

    template <class T, class Allocator>
    void ConvertSparsetoDense(
        const Vector<T, VectSparse, Allocator>& V_sparse,
        Vector<T, VectFull, Allocator>& V_dense);
    template <class T, class Allocator>
    void ConvertDenseToSparse(
        const Vector<T, VectFull, Allocator> V_dense,
        Vector<T, VectSparse, Allocator>& V_sparse);

    template <class T>
    void Fill(T value, Matrix<T, Symmetric, RowSymSparse>& M);

    template <class T0,
              class T1, class Prop1, class Storage1, class Allocator1,
              class T2, class Prop2, class Storage2, class Allocator2>
    void AddMatrixPosition(T0 c, int pi, int pj,
                           const Matrix<T1, Prop1, Storage1, Allocator1>& B,
                           Matrix<T2, Prop2, Storage2, Allocator2>& A);
#ifndef SWIG
    template <class T, template <class U> class Allocator>
    void GetRowPointer(const Matrix<T, General, RowMajor, Allocator<T> >& M,
                       int i, Vector<T, VectFull, Allocator<T> >& V);
#endif

#ifdef VERDANDI_WITH_PETSC
    template <class T0, class Allocator0, class Model>
    void Reallocate(Matrix<T0, General, PETScMPIDense, Allocator0>& A, int i,
                    int j, const Model& model);
#endif

    template <class T0, class Prop0, class Storage0, class Allocator0,
              class Model>
    void Reallocate(Matrix<T0, Prop0, Storage0, Allocator0>& A, int i, int j,
                    const Model& model);

#ifdef VERDANDI_WITH_PETSC
    template <class T0, class Allocator0, class Model>
    void Reallocate(Vector<T0, PETScPar, Allocator0>& V, int i,
                    const Model& model);
#endif

    template <class T0, class Storage0, class Allocator0, class Model>
    void Reallocate(Vector<T0, Storage0, Allocator0>& V, int i,
                    const Model& model);


} // namespace Verdandi.


#define VERDANDI_FILE_SHARE_USEFULFUNCTION_HXX
#endif
