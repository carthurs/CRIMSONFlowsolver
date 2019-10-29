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


#ifndef SELDON_FILE_MATRIX_CONVERSIONS_CXX


#include "Matrix_Conversions.hxx"


namespace Seldon
{

  /*
    From CSR formats to "Matlab" coordinate format.
    sym = true => the upper part and lower part are both generated
  */


  //! Conversion from RowSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    for (i = 0; i < m; i++)
      for (j = ptr[i]; j< ptr[i+1]; j++)
	{
	  IndRow(j) = i + index;
	  IndCol(j) = ind[j] + index;
	  Val(j) = val[j];
	}
  }


  //! Conversion from ColSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int n = A.GetN();
    int nnz = A.GetDataSize();
    IndCol.Reallocate(nnz);
    IndRow.Reallocate(nnz);
    Val.Reallocate(nnz);
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    for (i = 0; i < n; i++)
      for (j = ptr[i]; j< ptr[i+1]; j++)
	{
	  IndCol(j) = i + index;
	  IndRow(j) = ind[j] + index;
	  Val(j) = val[j];
	}
  }


  //! Conversion from RowSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    if (sym)
      {
	nnz *= 2;
	for (i = 0; i < m; i++)
	  if (ind[ptr[i]] == i)
	    nnz--;

	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = ind[j] + index;
	      Val(nb) = val[j];
	      Ptr(ind[j])++;
	      nb++;

	      if (ind[j] != i)
		{
		  IndRow(nb) = ind[j] + index;
		  IndCol(nb) = i + index;
		  Val(nb) = val[j];
		  Ptr(i)++;
		  nb++;
		}
	    }

	// Sorting the row numbers...
	Sort(IndRow, IndCol, Val);

	// ... and the column numbers.
	int offset = 0;
	for (i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }

      }
    else
      {
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j< ptr[i + 1]; j++)
	    {
	      IndRow(j) = i + index;
	      IndCol(j) = ind[j] + index;
	      Val(j) = val[j];
	    }
      }
  }


  //! Conversion from ColSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    if (sym)
      {
	nnz *= 2;
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    if (ind[j] == i)
	      nnz--;

	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = ind[j] + index;
	      Val(nb) = val[j];
	      Ptr(ind[j])++;
	      nb++;

	      if (ind[j] != i)
		{
		  IndRow(nb) = ind[j] + index;
		  IndCol(nb) = i + index;
		  Val(nb) = val[j];
		  Ptr(i)++;
		  nb++;
		}
	    }

	// Sorting the row numbers...
	Sort(IndRow, IndCol, Val);

	// ...and the column numbers.
	int offset = 0;
	for (i = 0; i < m; i++)
	  {
 	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }

      }
    else
      {
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j< ptr[i + 1]; j++)
	    {
	      IndRow(j) = i + index;
	      IndCol(j) = ind[j] + index;
	      Val(j) = val[j];
	    }
      }
  }


  //! Conversion from RowComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int m = A.GetM();
    int nnz = A.GetRealDataSize() + A.GetImagDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    nnz = 0;
    int* real_ptr = A.GetRealPtr();
    int* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    T* real_data = A.GetRealData();
    T* imag_data = A.GetImagData();
    IVect col; Vector<complex<T> > value;
    for (int i = 0; i < m; i++)
      {
        int nb_r = real_ptr[i+1] - real_ptr[i];
        int nb_i = imag_ptr[i+1] - imag_ptr[i];
        int size_row = nb_r + nb_i;
        if (size_row > col.GetM())
          {
            col.Reallocate(size_row);
            value.Reallocate(size_row);
          }

        int nb = 0;
        for (int j = real_ptr[i]; j < real_ptr[i+1]; j++)
          {
            col(nb) = real_ind[j] + index;
            value(nb) = complex<T>(real_data[j], 0);
            nb++;
          }

        for (int j = imag_ptr[i]; j < imag_ptr[i+1]; j++)
          {
            col(nb) = imag_ind[j] + index;
            value(nb) = complex<T>(0, imag_data[j]);
            nb++;
          }

        Assemble(nb, col, value);
        for (int j = 0; j < nb; j++)
          {
            IndRow(nnz + j) = index + i;
            IndCol(nnz + j) = col(j);
            Val(nnz + j) = value(j);
          }

        nnz += nb;
      }

    IndRow.Resize(nnz);
    IndCol.Resize(nnz);
    Val.Resize(nnz);
  }


  //! Conversion from ColComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int n = A.GetN();
    int nnz = A.GetRealDataSize() + A.GetImagDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    nnz = 0;
    int* real_ptr = A.GetRealPtr();
    int* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    T* real_data = A.GetRealData();
    T* imag_data = A.GetImagData();
    IVect col; Vector<complex<T> > value;
    for (int i = 0; i < n; i++)
      {
        int nb_r = real_ptr[i+1] - real_ptr[i];
        int nb_i = imag_ptr[i+1] - imag_ptr[i];
        int size_col = nb_r + nb_i;
        if (size_col > col.GetM())
          {
            col.Reallocate(size_col);
            value.Reallocate(size_col);
          }

        int nb = 0;
        for (int j = real_ptr[i]; j < real_ptr[i+1]; j++)
          {
            col(nb) = real_ind[j] + index;
            value(nb) = complex<T>(real_data[j], 0);
            nb++;
          }

        for (int j = imag_ptr[i]; j < imag_ptr[i+1]; j++)
          {
            col(nb) = imag_ind[j] + index;
            value(nb) = complex<T>(0, imag_data[j]);
            nb++;
          }

        Assemble(nb, col, value);
        for (int j = 0; j < nb; j++)
          {
            IndRow(nnz + j) = col(j);
            IndCol(nnz + j) = index + i;
            Val(nnz + j) = value(j);
          }

        nnz += nb;
      }

    IndRow.Resize(nnz);
    IndCol.Resize(nnz);
    Val.Resize(nnz);
  }


    //! Conversion from ArrayRowSymComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int m = A.GetM();
    int nnz = A.GetDataSize();
    int* real_ptr = A.GetRealPtr();
    int* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    T* real_data = A.GetRealData();
    T* imag_data = A.GetImagData();

    if (sym)
      {
	nnz *= 2;
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
        nnz = 0;
        IVect col; Vector<complex<T> > value;
	for (int i = 0; i < m; i++)
	  {
            int nb_r = real_ptr[i+1] - real_ptr[i];
            int nb_i = imag_ptr[i+1] - imag_ptr[i];
            int size_row = nb_r + nb_i;
            if (size_row > col.GetM())
              {
                col.Reallocate(size_row);
                value.Reallocate(size_row);
              }

            int nb = 0;
            for (int j = real_ptr[i]; j < real_ptr[i+1]; j++)
              {
                col(nb) = real_ind[j];
                value(nb) = complex<T>(real_data[j], 0);
                nb++;
              }

            for (int j = imag_ptr[i]; j < imag_ptr[i+1]; j++)
              {
                col(nb) = imag_ind[j];
                value(nb) = complex<T>(0, imag_data[j]);
                nb++;
              }

            Assemble(nb, col, value);
            for (int j = 0; j < nb; j++)
              {
                IndRow(nnz) = i + index;
                IndCol(nnz) = col(j) + index;
                Val(nnz) = value(j);
                Ptr(i)++;
                nnz++;

                if (col(j) != i)
                  {
                    IndRow(nnz) = col(j) + index;
                    IndCol(nnz) = i + index;
                    Val(nnz) = value(j);
                    Ptr(col(j))++;
                    nnz++;
                  }
              }
          }

        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);

        // Sorting the row numbers...
	Sort(IndRow, IndCol, Val);

	// ...and the column numbers.
	int offset = 0;
	for (int i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
        nnz = 0;
        IVect col; Vector<complex<T> > value;
	for (int i = 0; i < m; i++)
	  {
            int nb_r = real_ptr[i+1] - real_ptr[i];
            int nb_i = imag_ptr[i+1] - imag_ptr[i];
            int size_row = nb_r + nb_i;
            if (size_row > col.GetM())
              {
                col.Reallocate(size_row);
                value.Reallocate(size_row);
              }

            int nb = 0;
            for (int j = real_ptr[i]; j < real_ptr[i+1]; j++)
              {
                col(nb) = real_ind[j] + index;
                value(nb) = complex<T>(real_data[j], 0);
                nb++;
              }

            for (int j = real_ptr[i]; j < real_ptr[i+1]; j++)
              {
                col(nb) = imag_ind[j] + index;
                value(nb) = complex<T>(0, imag_ind[j]);
                nb++;
              }

            Assemble(nb, col, value);
            for (int j = 0; j < nb; j++)
              {
                IndRow(nnz + j) = index + i;
                IndCol(nnz + j) = col(j);
                Val(nnz + j) = value(j);
              }

            nnz += nb;
          }

        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);
      }
  }


  //! Conversion from ColSymComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int m = A.GetM();
    int nnz = A.GetDataSize();
    int* real_ptr = A.GetRealPtr();
    int* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    T* real_data = A.GetRealData();
    T* imag_data = A.GetImagData();

    if (sym)
      {
	nnz *= 2;
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
        nnz = 0;
        IVect row; Vector<complex<T> > value;
	for (int i = 0; i < m; i++)
	  {
            int nb_r = real_ptr[i+1] - real_ptr[i];
            int nb_i = imag_ptr[i+1] - imag_ptr[i];
            int size_col = nb_r + nb_i;
            if (size_col > row.GetM())
              {
                row.Reallocate(size_col);
                value.Reallocate(size_col);
              }

            int nb = 0;
            for (int j = real_ptr[i]; j < real_ptr[i+1]; j++)
              {
                row(nb) = real_ind[j];
                value(nb) = complex<T>(real_data[j], 0);
                nb++;
              }

            for (int j = imag_ptr[i]; j < imag_ptr[i+1]; j++)
              {
                row(nb) = imag_ind[j];
                value(nb) = complex<T>(0, imag_data[j]);
                nb++;
              }

            Assemble(nb, row, value);
            for (int j = 0; j < nb; j++)
              {
                IndRow(nnz) = i + index;
                IndCol(nnz) = row(j) + index;
                Val(nnz) = value(j);
                Ptr(i)++;
                nnz++;

                if (row(j) != i)
                  {
                    IndRow(nnz) = row(j) + index;
                    IndCol(nnz) = i + index;
                    Val(nnz) = value(j);
                    Ptr(row(j))++;
                    nnz++;
                  }
              }
          }

        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);

        // Sorting the column numbers...
	Sort(IndCol, IndRow, Val);

	// ...and the row numbers.
	int offset = 0;
	for (int i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndRow, Val);
	    offset += Ptr(i);
	  }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
        nnz = 0;
        IVect row; Vector<complex<T> > value;
	for (int i = 0; i < m; i++)
	  {
            int nb_r = real_ptr[i+1] - real_ptr[i];
            int nb_i = imag_ptr[i+1] - imag_ptr[i];
            int size_col = nb_r + nb_i;
            if (size_col > row.GetM())
              {
                row.Reallocate(size_col);
                value.Reallocate(size_col);
              }

            int nb = 0;
            for (int j = real_ptr[i]; j < real_ptr[i+1]; j++)
              {
                row(nb) = real_ind[j] + index;
                value(nb) = complex<T>(real_data[j], 0);
                nb++;
              }

            for (int j = imag_ptr[i]; j < imag_ptr[i+1]; j++)
              {
                row(nb) = imag_ind[j] + index;
                value(nb) = complex<T>(0, imag_data[j]);
                nb++;
              }

            Assemble(nb, row, value);
            for (int j = 0; j < nb; j++)
              {
                IndCol(nnz + j) = index + i;
                IndRow(nnz + j) = row(j);
                Val(nnz + j) = value(j);
              }

            nnz += nb;
          }

        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);
      }
  }


  /*
    From Sparse Array formats to "Matlab" coordinate format.
  */


  //! Conversion from ArrayRowSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    int nb = 0;
    for (i = 0; i < m; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	{
	  IndRow(nb) = i + index;
	  IndCol(nb) = A.Index(i, j) + index;
	  Val(nb) = A.Value(i, j);
	  nb++;
	}
  }


  //! Conversion from ArrayColSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int n = A.GetN();
    int nnz = A.GetDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    int nb = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetColumnSize(i); j++)
	{
	  IndRow(nb) = A.Index(i, j) + index;
	  IndCol(nb) = i + index;
	  Val(nb) = A.Value(i, j);
	  nb++;
	}
  }


  //! Conversion from ArrayRowComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int m = A.GetM();
    int nnz = A.GetRealDataSize() + A.GetImagDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    nnz = 0;
    IVect col; Vector<complex<T> > value;
    for (int i = 0; i < m; i++)
      {
        int size_row = A.GetRealRowSize(i) + A.GetImagRowSize(i);
        if (size_row > col.GetM())
          {
            col.Reallocate(size_row);
            value.Reallocate(size_row);
          }

        int nb = 0;
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            col(nb) = A.IndexReal(i, j) + index;
            value(nb) = complex<T>(A.ValueReal(i, j), 0);
            nb++;
          }

        for (int j = 0; j < A.GetImagRowSize(i); j++)
          {
            col(nb) = A.IndexImag(i, j) + index;
            value(nb) = complex<T>(0, A.ValueImag(i, j));
            nb++;
          }

        Assemble(nb, col, value);
        for (int j = 0; j < nb; j++)
          {
            IndRow(nnz + j) = index + i;
            IndCol(nnz + j) = col(j);
            Val(nnz + j) = value(j);
          }

        nnz += nb;
      }

    IndRow.Resize(nnz);
    IndCol.Resize(nnz);
    Val.Resize(nnz);
  }


  //! Conversion from ArrayColComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int nnz = A.GetRealDataSize() + A.GetImagDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    nnz = 0;
    IVect row; Vector<complex<T> > value;
    for (int i = 0; i < A.GetN(); i++)
      {
        int size_col = A.GetRealColumnSize(i) + A.GetImagColumnSize(i);
        if (size_col > row.GetM())
          {
            row.Reallocate(size_col);
            value.Reallocate(size_col);
          }

        int nb = 0;
        for (int j = 0; j < A.GetRealColumnSize(i); j++)
          {
            row(nb) = A.IndexReal(i, j) + index;
            value(nb) = complex<T>(A.ValueReal(i, j), 0);
            nb++;
          }

        for (int j = 0; j < A.GetImagColumnSize(i); j++)
          {
            row(nb) = A.IndexImag(i, j) + index;
            value(nb) = complex<T>(0, A.ValueImag(i, j));
            nb++;
          }

        Assemble(nb, row, value);
        for (int j = 0; j < nb; j++)
          {
            IndRow(nnz + j) = row(j);
            IndCol(nnz + j) = index + i;
            Val(nnz + j) = value(j);
          }

        nnz += nb;
      }

    IndRow.Resize(nnz);
    IndCol.Resize(nnz);
    Val.Resize(nnz);
  }


  //! Conversion from ArrayRowSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    if (sym)
      {
	nnz *= 2;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    if (A.Index(i, j) == i)
	      nnz--;

	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = A.Index(i, j) + index;
	      Val(nb) = A.Value(i, j);
	      Ptr(A.Index(i, j))++;
	      nb++;

	      if (A.Index(i, j) != i)
		{
		  IndRow(nb) = A.Index(i, j) + index;
		  IndCol(nb) = i + index;
		  Val(nb) = A.Value(i, j);
		  Ptr(i)++;
		  nb++;
		}
	    }

        // Sorting the row numbers...
	Sort(IndRow, IndCol, Val);

	// ...and the column numbers.
	int offset = 0;
	for (i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = A.Index(i, j) + index;
	      Val(nb) = A.Value(i, j);
	      nb++;
	    }
      }
  }


  //! Conversion from ArrayColSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int m = A.GetM();
    int nnz = A.GetDataSize();
    if (sym)
      {
	nnz *= 2;
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < A.GetColumnSize(i); j++)
	    if (A.Index(i, j) == i)
	      nnz--;

	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < A.GetColumnSize(i); j++)
	    {
	      IndCol(nb) = i + index;
	      IndRow(nb) = A.Index(i, j) + index;
	      Val(nb) = A.Value(i, j);
	      Ptr(A.Index(i, j))++;
	      nb++;

	      if (A.Index(i, j) != i)
		{
		  IndCol(nb) = A.Index(i, j) + index;
		  IndRow(nb) = i + index;
		  Val(nb) = A.Value(i, j);
		  Ptr(i)++;
		  nb++;
		}
	    }

	// Sorting the row numbers...
	Sort(IndRow, IndCol, Val);

	// ...and the column numbers.
	int offset = 0;
	for (int i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	int nb = 0;
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < A.GetColumnSize(i); j++)
	    {
	      IndRow(nb) = A.Index(i, j) + index;
	      IndCol(nb) = i + index;
	      Val(nb) = A.Value(i, j);
	      nb++;
	    }
      }
  }


  //! Conversion from ArrayRowSymComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int m = A.GetM();
    int nnz = A.GetDataSize();
    if (sym)
      {
	nnz *= 2;
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
        nnz = 0;
        IVect col; Vector<complex<T> > value;
	for (int i = 0; i < m; i++)
	  {
            int size_row = A.GetRealRowSize(i) + A.GetImagRowSize(i);
            if (size_row > col.GetM())
              {
                col.Reallocate(size_row);
                value.Reallocate(size_row);
              }

            int nb = 0;
            for (int j = 0; j < A.GetRealRowSize(i); j++)
              {
                col(nb) = A.IndexReal(i, j);
                value(nb) = complex<T>(A.ValueReal(i, j), 0);
                nb++;
              }

            for (int j = 0; j < A.GetImagRowSize(i); j++)
              {
                col(nb) = A.IndexImag(i, j);
                value(nb) = complex<T>(0, A.ValueImag(i, j));
                nb++;
              }

            Assemble(nb, col, value);
            for (int j = 0; j < nb; j++)
              {
                IndRow(nnz) = i + index;
                IndCol(nnz) = col(j) + index;
                Val(nnz) = value(j);
                Ptr(i)++;
                nnz++;

                if (col(j) != i)
                  {
                    IndRow(nnz) = col(j) + index;
                    IndCol(nnz) = i + index;
                    Val(nnz) = value(j);
                    Ptr(col(j))++;
                    nnz++;
                  }
              }
          }

        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);

        // Sorting the row numbers...
	Sort(IndRow, IndCol, Val);

	// ...and the column numbers.
	int offset = 0;
	for (int i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
        nnz = 0;
        IVect col; Vector<complex<T> > value;
	for (int i = 0; i < m; i++)
	  {
            int size_row = A.GetRealRowSize(i) + A.GetImagRowSize(i);
            if (size_row > col.GetM())
              {
                col.Reallocate(size_row);
                value.Reallocate(size_row);
              }

            int nb = 0;
            for (int j = 0; j < A.GetRealRowSize(i); j++)
              {
                col(nb) = A.IndexReal(i, j) + index;
                value(nb) = complex<T>(A.ValueReal(i, j), 0);
                nb++;
              }

            for (int j = 0; j < A.GetImagRowSize(i); j++)
              {
                col(nb) = A.IndexImag(i, j) + index;
                value(nb) = complex<T>(0, A.ValueImag(i, j));
                nb++;
              }

            Assemble(nb, col, value);
            for (int j = 0; j < nb; j++)
              {
                IndRow(nnz + j) = index + i;
                IndCol(nnz + j) = col(j);
                Val(nnz + j) = value(j);
              }

            nnz += nb;
          }

        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);
      }
  }


  //! Conversion from ArrayColSymComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int m = A.GetM();
    int nnz = A.GetDataSize();
    if (sym)
      {
	nnz *= 2;
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
        nnz = 0;
        IVect row; Vector<complex<T> > value;
	for (int i = 0; i < m; i++)
	  {
            int size_col = A.GetRealColumnSize(i) + A.GetImagColumnSize(i);
            if (size_col > row.GetM())
              {
                row.Reallocate(size_col);
                value.Reallocate(size_col);
              }

            int nb = 0;
            for (int j = 0; j < A.GetRealColumnSize(i); j++)
              {
                row(nb) = A.IndexReal(i, j);
                value(nb) = complex<T>(A.ValueReal(i, j), 0);
                nb++;
              }

            for (int j = 0; j < A.GetImagColumnSize(i); j++)
              {
                row(nb) = A.IndexImag(i, j);
                value(nb) = complex<T>(0, A.ValueImag(i, j));
                nb++;
              }

            Assemble(nb, row, value);
            for (int j = 0; j < nb; j++)
              {
                IndRow(nnz) = i + index;
                IndCol(nnz) = row(j) + index;
                Val(nnz) = value(j);
                Ptr(i)++;
                nnz++;

                if (row(j) != i)
                  {
                    IndRow(nnz) = row(j) + index;
                    IndCol(nnz) = i + index;
                    Val(nnz) = value(j);
                    Ptr(row(j))++;
                    nnz++;
                  }
              }
          }

        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);

        // Sorting the column numbers...
	Sort(IndCol, IndRow, Val);

	// ...and the row numbers.
	int offset = 0;
	for (int i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndRow, Val);
	    offset += Ptr(i);
	  }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
        nnz = 0;
        IVect row; Vector<complex<T> > value;
	for (int i = 0; i < m; i++)
	  {
            int size_col = A.GetRealColumnSize(i) + A.GetImagColumnSize(i);
            if (size_col > row.GetM())
              {
                row.Reallocate(size_col);
                value.Reallocate(size_col);
              }

            int nb = 0;
            for (int j = 0; j < A.GetRealColumnSize(i); j++)
              {
                row(nb) = A.IndexReal(i, j) + index;
                value(nb) = complex<T>(A.ValueReal(i, j), 0);
                nb++;
              }

            for (int j = 0; j < A.GetImagColumnSize(i); j++)
              {
                row(nb) = A.IndexImag(i, j) + index;
                value(nb) = complex<T>(0, A.ValueImag(i, j));
                nb++;
              }

            Assemble(nb, row, value);
            for (int j = 0; j < nb; j++)
              {
                IndCol(nnz + j) = index + i;
                IndRow(nnz + j) = row(j);
                Val(nnz + j) = value(j);
              }

            nnz += nb;
          }

        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);
      }
  }


  /*
    From "Matlab" coordinate format to CSR formats.
  */


  //! Conversion from coordinate format to RowSparse.
  /*! Contrary to the other conversion functions
    ConvertMatrix_from_Coordinates, this one accepts duplicates.
    \param[in] IndRow_ row indexes of the non-zero elements.
    \param[in] IndCol_ column indexes of the non-zero elements.
    \param[in] Val values of the non-zero elements.
    \param[out] A matrix defined by \a IndRow, \a IndCol and \a Val.
    \param[in] index index of the first column and of the first row.
  */
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowSparse, Allocator3>& A,
				 int index)
  {
    int Nelement = IndRow_.GetLength();

    Vector<int, VectFull, CallocAlloc<int> > IndRow(Nelement),
      IndCol(Nelement);

    for (int i = 0; i < Nelement; i++)
      {
	IndRow(i) = IndRow_(i) - index;
	IndCol(i) = IndCol_(i) - index;
      }

    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();

    int m = row_max + 1;
    int n = col_max + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());

    Sort(IndRow, IndCol, Val);

    // Construction of array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(m + 1);
    Ptr.Zero();
    for (int i = 0; i < Nelement; i++)
      Ptr(IndRow(i)+1)++;

    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);

    // Sorts 'IndCol'
    for (int i = 0; i < m; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndCol, Val);

    A.SetData(m, n, Val, Ptr, IndCol);
  }


  //! Conversion from coordinate format to ColSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColSparse, Allocator3>& A,
				 int index)
  {
    // Assuming that there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());

    // Sorts the array 'IndCol'.
    Sort(IndCol, IndRow, Val);

    // Construction of array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(n + 1);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndCol(i) + 1)++;
      }

    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);

    // Sorts 'IndRow'
    for (int i = 0; i < n; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndRow, Val);

    A.SetData(m, n, Val, Ptr, IndRow);
  }


  //! Conversion from coordinate format to RowSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowSymSparse, Allocator3>& A,
				 int index)
  {
    // Assuming there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (int i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;

    if (nb_low > 0)
      {
	int nb = 0;
	for (int i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }

	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Construction of array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(m + 1);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i) + 1)++;
      }

    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);

    // Sorts 'IndCol'.
    for (int i = 0; i < m; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndCol, Val);

    A.SetData(m, n, Val, Ptr, IndCol);
  }


  //! Conversion from coordinate format to ColSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColSymSparse, Allocator3>& A,
				 int index)
  {
    // Assuming there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }

    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (int i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;

    if (nb_low > 0)
      {
	int nb = 0;
	for (int i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }

	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndCol'.
    Sort(IndCol, IndRow, Val);

    // Construction of array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(m + 1);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndCol(i) + 1)++;
      }

    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);

    // Sorts 'IndRow'.
    for (int i = 0; i < m; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndRow, Val);

    A.SetData(m, n, Val, Ptr, IndRow);
  }


  //! Conversion from coordinate format to RowComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowComplexSparse,
				 Allocator4>& A,
				 int index)
  {
    if (IndRow.GetM() <= 0)
      {
        A.Clear();
        return;
      }

    T zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Number of elements per row.
    Vector<int, VectFull, CallocAlloc<int> > PtrReal(m+1), PtrImag(m+1), Ptr(m);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndRow.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndRow(i))++;
	if (real(Val(i)) != zero)
          PtrReal(IndRow(i)+1)++;

	if (imag(Val(i)) != zero)
          PtrImag(IndRow(i)+1)++;
      }

    for (int i = 0; i < m; i++)
      {
        PtrReal(i+1) += PtrReal(i);
        PtrImag(i+1) += PtrImag(i);
      }
    int real_nz = PtrReal(m), imag_nz = PtrImag(m);

    // Fills matrix 'A'.
    Vector<int, VectFull, CallocAlloc<int> > IndReal(real_nz), IndImag(imag_nz);
    Vector<T, VectFull, Allocator4> ValReal(real_nz), ValImag(imag_nz);
    int offset = 0;
    for (int i = 0; i < m; i++)
      {
        int nb = PtrReal(i);
        for (int j = 0; j < Ptr(i); j++)
          if (real(Val(offset + j)) != zero)
            {
              IndReal(nb) = IndCol(offset + j);
              ValReal(nb) = real(Val(offset + j));
              nb++;
            }

        nb = PtrImag(i);
        for (int j = 0; j < Ptr(i); j++)
          if (imag(Val(offset + j)) != zero)
            {
              IndImag(nb) = IndCol(offset + j);
              ValImag(nb) = imag(Val(offset + j));
              nb++;
            }

        // sorting column numbers
        Sort(PtrReal(i), PtrReal(i+1)-1, IndReal, ValReal);
        Sort(PtrImag(i), PtrImag(i+1)-1, IndImag, ValImag);

        offset += Ptr(i);
      }

    // providing pointers to A
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }


  //! Conversion from coordinate format to ColComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColComplexSparse,
				 Allocator4>& A,
				 int index)
  {
    if (IndRow.GetM() <= 0)
      {
        A.Clear();
        return;
      }

    T zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // Sorts the array 'IndCol'.
    Sort(IndCol, IndRow, Val);

    // Number of elements per column.
    Vector<int, VectFull, CallocAlloc<int> > PtrReal(n+1), PtrImag(n+1), Ptr(n);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndCol.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndCol(i))++;
	if (real(Val(i)) != zero)
          PtrReal(IndCol(i)+1)++;

	if (imag(Val(i)) != zero)
          PtrImag(IndCol(i)+1)++;
      }

    for (int i = 0; i < n; i++)
      {
        PtrReal(i+1) += PtrReal(i);
        PtrImag(i+1) += PtrImag(i);
      }
    int real_nz = PtrReal(n), imag_nz = PtrImag(n);

    // Fills matrix 'A'.
    Vector<int, VectFull, CallocAlloc<int> > IndReal(real_nz), IndImag(imag_nz);
    Vector<T, VectFull, Allocator4> ValReal(real_nz), ValImag(imag_nz);
    int offset = 0;
    for (int i = 0; i < n; i++)
      {
        int nb = PtrReal(i);
        for (int j = 0; j < Ptr(i); j++)
          if (real(Val(offset + j)) != zero)
            {
              IndReal(nb) = IndRow(offset + j);
              ValReal(nb) = real(Val(offset + j));
              nb++;
            }

        nb = PtrImag(i);
        for (int j = 0; j < Ptr(i); j++)
          if (imag(Val(offset + j)) != zero)
            {
              IndImag(nb) = IndRow(offset + j);
              ValImag(nb) = imag(Val(offset + j));
              nb++;
            }

        // sorting column numbers
        Sort(PtrReal(i), PtrReal(i+1)-1, IndReal, ValReal);
        Sort(PtrImag(i), PtrImag(i+1)-1, IndImag, ValImag);

        offset += Ptr(i);
      }

    // providing pointers to A
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }


  //! Conversion from coordinate format to RowSymComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowSymComplexSparse,
				 Allocator4>& A,
				 int index)
  {
    if (IndRow.GetM() <= 0)
      {
        A.Clear();
        return;
      }

    T zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Number of elements per row.
    Vector<int, VectFull, CallocAlloc<int> > PtrReal(m+1), PtrImag(m+1), Ptr(m);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndRow.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndRow(i))++;
	if (IndRow(i) <= IndCol(i))
          {
            if (real(Val(i)) != zero)
              PtrReal(IndRow(i)+1)++;

            if (imag(Val(i)) != zero)
              PtrImag(IndRow(i)+1)++;
          }
      }

    for (int i = 0; i < m; i++)
      {
        PtrReal(i+1) += PtrReal(i);
        PtrImag(i+1) += PtrImag(i);
      }

    int real_nz = PtrReal(m), imag_nz = PtrImag(m);

    // Fills matrix 'A'.
    Vector<int, VectFull, CallocAlloc<int> > IndReal(real_nz), IndImag(imag_nz);
    Vector<T, VectFull, Allocator4> ValReal(real_nz), ValImag(imag_nz);
    int offset = 0;
    for (int i = 0; i < m; i++)
      {
        int nb = PtrReal(i);
        for (int j = 0; j < Ptr(i); j++)
          if (i <= IndCol(offset+j))
            if (real(Val(offset + j)) != zero)
              {
                IndReal(nb) = IndCol(offset + j);
                ValReal(nb) = real(Val(offset + j));
                nb++;
              }

        nb = PtrImag(i);
        for (int j = 0; j < Ptr(i); j++)
          if (i <= IndCol(offset+j))
            if (imag(Val(offset + j)) != zero)
              {
                IndImag(nb) = IndCol(offset + j);
                ValImag(nb) = imag(Val(offset + j));
                nb++;
            }

        // sorting column numbers
        Sort(PtrReal(i), PtrReal(i+1)-1, IndReal, ValReal);
        Sort(PtrImag(i), PtrImag(i+1)-1, IndImag, ValImag);

        offset += Ptr(i);
      }

    // providing pointers to A
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }


  //! Conversion from coordinate format to ColSymComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColSymComplexSparse,
				 Allocator4>& A,
				 int index)
  {
    if (IndRow.GetM() <= 0)
      {
        A.Clear();
        return;
      }

    T zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // Sorts the array 'IndCol'.
    Sort(IndCol, IndRow, Val);

    // Number of elements per column.
    Vector<int, VectFull, CallocAlloc<int> > PtrReal(n+1), PtrImag(n+1), Ptr(n);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndCol.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndCol(i))++;
	if (IndRow(i) <= IndCol(i))
          {
            if (real(Val(i)) != zero)
              PtrReal(IndCol(i)+1)++;

            if (imag(Val(i)) != zero)
              PtrImag(IndCol(i)+1)++;
          }
      }

    for (int i = 0; i < n; i++)
      {
        PtrReal(i+1) += PtrReal(i);
        PtrImag(i+1) += PtrImag(i);
      }

    int real_nz = PtrReal(n), imag_nz = PtrImag(n);

    // Fills matrix 'A'.
    Vector<int, VectFull, CallocAlloc<int> > IndReal(real_nz), IndImag(imag_nz);
    Vector<T, VectFull, Allocator4> ValReal(real_nz), ValImag(imag_nz);
    int offset = 0;
    for (int i = 0; i < n; i++)
      {
        int nb = PtrReal(i);
        for (int j = 0; j < Ptr(i); j++)
          if (IndRow(offset+j) <= i)
            if (real(Val(offset + j)) != zero)
              {
                IndReal(nb) = IndRow(offset + j);
                ValReal(nb) = real(Val(offset + j));
                nb++;
              }

        nb = PtrImag(i);
        for (int j = 0; j < Ptr(i); j++)
          if (IndRow(offset+j) <= i)
            if (imag(Val(offset + j)) != zero)
              {
                IndImag(nb) = IndRow(offset + j);
                ValImag(nb) = imag(Val(offset + j));
                nb++;
              }

        // sorting column numbers
        Sort(PtrReal(i), PtrReal(i+1)-1, IndReal, ValReal);
        Sort(PtrImag(i), PtrImag(i+1)-1, IndImag, ValImag);

        offset += Ptr(i);
      }

    // providing pointers to A
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }


  /*
    From Sparse Array formats to "Matlab" coordinate format.
  */


  //! Conversion from coordinate format to ArrayRowSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowSparse,
				 Allocator3>& A,
				 int index)
  {
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Number of elements per row.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(m);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i))++;
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      if (Ptr(i) > 0)
	{
	  A.ReallocateRow(i, Ptr(i));
	  for (int j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndCol(offset + j);
	      A.Value(i, j) = Val(offset + j);
	    }
	  offset += Ptr(i);
	}

    // Assembles 'A' to sort column numbers.
    A.Assemble();
  }


  //! Conversion from coordinate format to ArrayColSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColSparse,
				 Allocator3>& A,
				 int index)
  {
    // Assuming there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // Sorts array 'IndCol'.
    Sort(IndCol, IndRow, Val);

    // Construction of array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(n);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndCol(i))++;
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < n; i++)
      if (Ptr(i) > 0)
	{
	  A.ReallocateColumn(i, Ptr(i));
	  for (int j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndRow(offset + j);
	      A.Value(i, j) = Val(offset + j);
	    }
	  offset += Ptr(i);
	}

    // Assembles 'A' to sort row numbers.
    A.Assemble();
  }


  //! Conversion from coordinate format to ArrayRowSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowSymSparse,
				 Allocator3>& A,
				 int index)
  {
    // Assuming that there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (int i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;

    if (nb_low > 0)
      {
	int nb = 0;
	for (int i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }

	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Construction of array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(m);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i))++;
      }

    // Fills matrix 'A'.
    A.Clear(); A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      if (Ptr(i) > 0)
	{
          // sorting column numbers
          Sort(offset, offset+Ptr(i)-1, IndCol, Val);

          // putting values in A
	  A.ReallocateRow(i, Ptr(i));
	  for (int j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndCol(offset + j);
	      A.Value(i, j) = Val(offset + j);
            }

	  offset += Ptr(i);
	}
  }


  //! Conversion from coordinate format to ArrayColSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColSymSparse,
				 Allocator3>& A,
				 int index)
  {
    // Assuming that there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }

    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (int i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;

    if (nb_low > 0)
      {
	int nb = 0;
	for (int i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }

	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndRow'.
    Sort(IndCol, IndRow, Val);

    // Construction of array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(m);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndCol(i))++;
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      if (Ptr(i) > 0)
	{
	  A.ReallocateColumn(i, Ptr(i));
	  for (int j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndRow(offset + j);
	      A.Value(i, j) = Val(offset + j);
	    }
	  offset += Ptr(i);
	}

    // Assembles 'A' to sort row numbers.
    A.Assemble();
  }


  //! Conversion from coordinate format to ArrayRowComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowComplexSparse,
				 Allocator4>& A,
				 int index)
  {
    if (IndRow.GetM() <= 0)
      {
        A.Clear();
        return;
      }

    T zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Number of elements per row.
    Vector<int> PtrReal(m), PtrImag(m), Ptr(m);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndRow.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndRow(i))++;
	if (real(Val(i)) != zero)
          PtrReal(IndRow(i))++;

	if (imag(Val(i)) != zero)
          PtrImag(IndRow(i))++;
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      {
        if (PtrReal(i) > 0)
          {
            A.ReallocateRealRow(i, PtrReal(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (real(Val(offset + j)) != zero)
                {
                  A.IndexReal(i, nb) = IndCol(offset + j);
                  A.ValueReal(i, nb) = real(Val(offset + j));
                  nb++;
                }
          }

        if (PtrImag(i) > 0)
          {
            A.ReallocateImagRow(i, PtrImag(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (imag(Val(offset + j)) != zero)
                {
                  A.IndexImag(i, nb) = IndCol(offset + j);
                  A.ValueImag(i, nb) = imag(Val(offset + j));
                  nb++;
                }
          }

        offset += Ptr(i);
      }

    // Assembles 'A' to sort column numbers.
    A.Assemble();
  }


  //! Conversion from coordinate format to ArrayColComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColComplexSparse,
				 Allocator4>& A,
				 int index)
  {
    if (IndRow.GetM() <= 0)
      {
        A.Clear();
        return;
      }

    T zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // Sorts the array 'IndRow'.
    Sort(IndCol, IndRow, Val);

    // Number of elements per row.
    Vector<int> PtrReal(n), PtrImag(n), Ptr(n);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndRow.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndCol(i))++;
	if (real(Val(i)) != zero)
          PtrReal(IndCol(i))++;

	if (imag(Val(i)) != zero)
          PtrImag(IndCol(i))++;
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < n; i++)
      {
        if (PtrReal(i) > 0)
          {
            A.ReallocateRealColumn(i, PtrReal(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (real(Val(offset + j)) != zero)
                {
                  A.IndexReal(i, nb) = IndRow(offset + j);
                  A.ValueReal(i, nb) = real(Val(offset + j));
                  nb++;
                }
          }

        if (PtrImag(i) > 0)
          {
            A.ReallocateImagColumn(i, PtrImag(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (imag(Val(offset + j)) != zero)
                {
                  A.IndexImag(i, nb) = IndRow(offset + j);
                  A.ValueImag(i, nb) = imag(Val(offset + j));
                  nb++;
                }
          }

        offset += Ptr(i);
      }

    // Assembles 'A' to sort row numbers.
    A.Assemble();
  }


  //! Conversion from coordinate format to ArrayRowSymComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowSymComplexSparse,
				 Allocator4>& A, int index)
  {
    if (IndRow.GetM() <= 0)
      {
        A.Clear();
        return;
      }

    T zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Number of elements per row.
    Vector<int> PtrReal(m), PtrImag(m), Ptr(m);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndRow.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndRow(i))++;
        if (IndRow(i) <= IndCol(i))
          {
            if (real(Val(i)) != zero)
              PtrReal(IndRow(i))++;

            if (imag(Val(i)) != zero)
              PtrImag(IndRow(i))++;
          }
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      {
        if (PtrReal(i) > 0)
          {
            A.ReallocateRealRow(i, PtrReal(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (real(Val(offset + j)) != zero)
                {
                  if (i <= IndCol(offset+j))
                    {
                      A.IndexReal(i, nb) = IndCol(offset + j);
                      A.ValueReal(i, nb) = real(Val(offset + j));
                      nb++;
                    }
                }
          }

        if (PtrImag(i) > 0)
          {
            A.ReallocateImagRow(i, PtrImag(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (imag(Val(offset + j)) != zero)
                {
                  if (i <= IndCol(offset+j))
                    {
                      A.IndexImag(i, nb) = IndCol(offset + j);
                      A.ValueImag(i, nb) = imag(Val(offset + j));
                      nb++;
                    }
                }
          }

        offset += Ptr(i);
      }

    // Assembles 'A' to sort column numbers.
    A.Assemble();
  }


  //! Conversion from coordinate format to ArrayRowSymComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColSymComplexSparse,
				 Allocator4>& A, int index)
  {
    if (IndRow.GetM() <= 0)
      {
        A.Clear();
        return;
      }

    T zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // Sorts the array 'IndRow'.
    Sort(IndCol, IndRow, Val);

    // Number of elements per row.
    Vector<int> PtrReal(m), PtrImag(m), Ptr(m);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndRow.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndCol(i))++;
        if (IndRow(i) <= IndCol(i))
          {
            if (real(Val(i)) != zero)
              PtrReal(IndCol(i))++;

            if (imag(Val(i)) != zero)
              PtrImag(IndCol(i))++;
          }
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      {
        if (PtrReal(i) > 0)
          {
            A.ReallocateRealColumn(i, PtrReal(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (real(Val(offset + j)) != zero)
                {
                  if (IndRow(offset+j) <= i)
                    {
                      A.IndexReal(i, nb) = IndRow(offset + j);
                      A.ValueReal(i, nb) = real(Val(offset + j));
                      nb++;
                    }
                }
          }

        if (PtrImag(i) > 0)
          {
            A.ReallocateImagColumn(i, PtrImag(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (imag(Val(offset + j)) != zero)
                {
                  if (IndRow(offset+j) <= i)
                    {
                      A.IndexImag(i, nb) = IndRow(offset + j);
                      A.ValueImag(i, nb) = imag(Val(offset + j));
                      nb++;
                    }
                }
          }

        offset += Ptr(i);
      }

    // Assembles 'A' to sort row numbers.
    A.Assemble();
  }


  /*
    From Sparse formats to CSC format
  */


  //! Conversion from RowSparse to CSC format
  /*!
    if sym_pat is equal to true, the pattern is symmetrized
    by adding artificial null entries
   */
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with nnz entries.
    int nnz = A.GetDataSize();
    int n = A.GetN();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();

    // Conversion in coordinate format.
    Vector<Tint, VectFull, CallocAlloc<Tint> > IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);

    // Sorting with respect to column numbers.
    Sort(IndCol, IndRow, Val);

    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(n + 1);
    Ptr.Fill(0);

    // Counting non-zero entries per column.
    for (int i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;

    int nb_new_val = 0;

    if (sym_pat)
      {
        // Counting entries that are on the symmetrized pattern without being
        // in the original pattern.
        int k = 0;
        for (int i = 0; i < n; i++)
          {
            while (k < IndCol.GetM() && IndCol(k) < i)
              k++;

            for (int j = ptr_[i]; j < ptr_[i+1]; j++)
              {
                int irow = ind_[j];
                while (k < IndCol.GetM() && IndCol(k) == i
                       && IndRow(k) < irow)
                  k++;

                if (k < IndCol.GetM() && IndCol(k) == i && IndRow(k) == irow)
                  // Already existing entry.
                  k++;
                else
                  {
                    // New entry.
                    Ptr(i + 1)++;
                    nb_new_val++;
                  }
              }
          }
      }

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);

    if (sym_pat && (nb_new_val > 0))
      {
        // Changing 'IndRow' and 'Val', and assembling the pattern.
        Vector<Tint, VectFull, Alloc3> OldInd(IndRow);
        Vector<T, VectFull, Alloc4> OldVal(Val);
        IndRow.Reallocate(nnz + nb_new_val);
        Val.Reallocate(nnz + nb_new_val);
        int k = 0, nb = 0;
        for (int i = 0; i < n; i++)
          {
            while (k < IndCol.GetM() && IndCol(k) < i)
              {
                IndRow(nb) = OldInd(k);
                Val(nb) = OldVal(k);
                nb++;
                k++;
              }

            for (int j = ptr_[i]; j < ptr_[i+1]; j++)
              {
                int irow = ind_[j];
                while (k < IndCol.GetM() && IndCol(k) == i
                       && OldInd(k) < irow)
                  {
                    IndRow(nb) = OldInd(k);
                    Val(nb) = OldVal(k);
                    nb++;
                    k++;
                  }

                if (k < IndCol.GetM() && IndCol(k) == i && OldInd(k) == irow)
                  {
                    // Already existing entry.
                    IndRow(nb) = OldInd(k);
                    Val(nb) = OldVal(k);
                    nb++;
                    k++;
                  }
                else
                  {
                    // New entry (null).
                    IndRow(nb) = irow;
                    Val(nb) = 0;
                    nb++;
                  }
              }
          }
      }
  }


  //! Conversion from ArrayRowSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with nnz entries.
    int nnz = A.GetDataSize();
    int n = A.GetN();

    // Conversion in coordinate format.
    Vector<Tint, VectFull, CallocAlloc<Tint> > IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);

    // Sorting with respect to column numbers.
    Sort(IndCol, IndRow, Val);

    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(n + 1);
    Ptr.Fill(0);

    // Counting non-zero entries per column.
    for (int i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;

    int nb_new_val = 0;

    if (sym_pat)
      {
        // Counting entries that are on the symmetrized pattern without being
        // in the original pattern.
        int k = 0;
        for (int i = 0; i < n; i++)
          {
            while (k < IndCol.GetM() && IndCol(k) < i)
              k++;

            for (int j = 0; j < A.GetRowSize(i); j++)
              {
                int irow = A.Index(i, j);
                while (k < IndCol.GetM() && IndCol(k) == i
                       && IndRow(k) < irow)
                  k++;

                if (k < IndCol.GetM() && IndCol(k) == i && IndRow(k) == irow)
                  // Already existing entry.
                  k++;
                else
                  {
                    // New entry.
                    Ptr(i + 1)++;
                    nb_new_val++;
                  }
              }
          }
      }

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);

    if (sym_pat && (nb_new_val > 0))
      {
        // Changing 'IndRow' and 'Val', and assembling the pattern.
        Vector<Tint, VectFull, Alloc3> OldInd(IndRow);
        Vector<T, VectFull, Alloc4> OldVal(Val);
        IndRow.Reallocate(nnz + nb_new_val);
        Val.Reallocate(nnz + nb_new_val);
        int k = 0, nb = 0;
        for (int i = 0; i < n; i++)
          {
            while (k < IndCol.GetM() && IndCol(k) < i)
              {
                IndRow(nb) = OldInd(k);
                Val(nb) = OldVal(k);
                nb++;
                k++;
              }

            for (int j = 0; j < A.GetRowSize(i); j++)
              {
                int irow = A.Index(i, j);
                while (k < IndCol.GetM() && IndCol(k) == i
                       && OldInd(k) < irow)
                  {
                    IndRow(nb) = OldInd(k);
                    Val(nb) = OldVal(k);
                    nb++;
                    k++;
                  }

                if (k < IndCol.GetM() && IndCol(k) == i && OldInd(k) == irow)
                  {
                    // Already existing entry.
                    IndRow(nb) = OldInd(k);
                    Val(nb) = OldVal(k);
                    nb++;
                    k++;
                  }
                else
                  {
                    // New entry (null).
                    IndRow(nb) = irow;
                    Val(nb) = 0;
                    nb++;
                  }
              }
          }
      }
  }


  //! Conversion from ColSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with 'nnz' entries.
    int nnz = A.GetDataSize();
    int n = A.GetN();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();

    // Conversion in coordinate format.
    Vector<Tint, VectFull, CallocAlloc<Tint> > IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);

    // Sorting with respect to row numbers.
    Sort(IndRow, IndCol, Val);

    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(n + 1);
    Ptr.Fill(0);

    // Counting non-zero entries per column.
    for (int i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;

    int nb_new_val = 0;

    if (sym_pat)
      {
        // Counting entries that are on the symmetrized pattern without being
        // in the original pattern.
        int k = 0;
        for (int i = 0; i < n; i++)
          {
            while (k < IndRow.GetM() && IndRow(k) < i)
              k++;

            for (int j = ptr_[i]; j < ptr_[i+1]; j++)
              {
                int icol = ind_[j];
                while (k < IndRow.GetM() && IndRow(k) == i
                       && IndCol(k) < icol)
                  k++;

                if (k < IndRow.GetM() && IndRow(k) == i && IndCol(k) == icol)
                  // Already existing entry.
                  k++;
                else
                  {
                    // New entry.
                    Ptr(i + 1)++;
                    nb_new_val++;
                  }
              }
          }
      }

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);

    if (sym_pat && (nb_new_val > 0))
      {
        // Changing 'IndRow' and 'Val', and assembling the pattern.
        Vector<Tint, VectFull, Alloc3> OldInd(IndRow);
        Vector<T, VectFull, Alloc4> OldVal(Val);
        IndRow.Reallocate(nnz + nb_new_val);
        Val.Reallocate(nnz + nb_new_val);
        int k = 0, nb = 0;
        for (int i = 0; i < n; i++)
          {
            while (k < OldInd.GetM() && OldInd(k) < i)
              {
		// null entries (due to symmetrisation)
                IndRow(nb) = IndCol(k);
                Val(nb) = 0;
                nb++;
                k++;
              }

            for (int j = ptr_[i]; j < ptr_[i+1]; j++)
              {
                int irow = ind_[j];
                while (k < OldInd.GetM() && OldInd(k) == i
                       && IndCol(k) < irow)
                  {
		    // null entries (due to symmetrisation)
                    IndRow(nb) = IndCol(k);
                    Val(nb) = 0;
                    nb++;
                    k++;
                  }

                if (k < OldInd.GetM() && OldInd(k) == i && IndCol(k) == irow)
                  {
                    // Already existing entry.
                    IndRow(nb) = IndCol(k);
                    Val(nb) = data_[j];
                    nb++;
                    k++;
                  }
                else
                  {
                    // New entry
                    IndRow(nb) = irow;
                    Val(nb) = data_[j];
                    nb++;
                  }
              }
          }
      }
    else
      {
	// sorting by columns
	Sort(IndCol, IndRow, Val);
      }

  }


  //! Conversion from ArrayColSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with 'nnz' entries.
    int nnz = A.GetDataSize();
    int n = A.GetN();

    // Conversion in coordinate format.
    Vector<Tint, VectFull, CallocAlloc<Tint> > IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);

    // Sorting with respect to row numbers.
    Sort(IndRow, IndCol, Val);

    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(n + 1);
    Ptr.Fill(0);

    // Counting non-zero entries per column.
    for (int i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;

    int nb_new_val = 0;

    if (sym_pat)
      {
        // Counting entries that are on the symmetrized pattern without being
        // in the original pattern.
        int k = 0;
        for (int i = 0; i < n; i++)
          {
            while (k < IndRow.GetM() && IndRow(k) < i)
              k++;

            for (int j = 0; j < A.GetColumnSize(i); j++)
              {
                int icol = A.Index(i, j);
                while (k < IndRow.GetM() && IndRow(k) == i
                       && IndCol(k) < icol)
                  k++;

                if (k < IndRow.GetM() && IndRow(k) == i && IndCol(k) == icol)
                  // Already existing entry.
                  k++;
                else
                  {
                    // New entry.
                    Ptr(i + 1)++;
                    nb_new_val++;
                  }
              }
          }
      }

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);

    if (sym_pat && (nb_new_val > 0))
      {
        // Changing 'IndRow' and 'Val', and assembling the pattern.
        Vector<Tint, VectFull, Alloc3> OldInd(IndRow);
        Vector<T, VectFull, Alloc4> OldVal(Val);
        IndRow.Reallocate(nnz + nb_new_val);
        Val.Reallocate(nnz + nb_new_val);
        int k = 0, nb = 0;
        for (int i = 0; i < n; i++)
          {
            while (k < OldInd.GetM() && OldInd(k) < i)
              {
		// null entries (due to symmetrisation)
                IndRow(nb) = IndCol(k);
                Val(nb) = 0;
                nb++;
                k++;
              }

            for (int j = 0; j < A.GetColumnSize(i); j++)
              {
                int irow = A.Index(i, j);
                while (k < OldInd.GetM() && OldInd(k) == i
                       && IndCol(k) < irow)
                  {
		    // null entries (due to symmetrisation)
                    IndRow(nb) = IndCol(k);
                    Val(nb) = 0;
                    nb++;
                    k++;
                  }

                if (k < OldInd.GetM() && OldInd(k) == i && IndCol(k) == irow)
                  {
                    // Already existing entry.
                    IndRow(nb) = IndCol(k);
                    Val(nb) = A.Value(i, j);
                    nb++;
                    k++;
                  }
                else
                  {
                    // New entry
                    IndRow(nb) = irow;
                    Val(nb) = A.Value(i, j);
                    nb++;
                  }
              }
          }
      }
    else
      {
	// sorting by columns
	Sort(IndCol, IndRow, Val);
      }

  }


  //! Conversion from ColSymSparse to symmetric CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& Ind,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();
    int nnz = A.GetDataSize();

    Ptr.Reallocate(n+1);
    Ind.Reallocate(nnz);
    Value.Reallocate(nnz);

    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();
    for (int i = 0; i <= n; i++)
      Ptr(i) = ptr_[i];

    for (int i = 0; i < nnz; i++)
      {
	Ind(i) = ind_[i];
	Value(i) = data_[i];
      }
  }


  //! Conversion from ColSymSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();

    Vector<Tint, VectFull, Alloc3> IndCol;

    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);

    // sorting by columns
    Sort(IndCol, IndRow, Value);

    Ptr.Reallocate(n+1);
    Ptr.Zero();
    // counting number of non-zero entries
    int nnz = 0;
    for (int i = 0; i < IndCol.GetM(); i++)
      {
	Ptr(IndCol(i) + 1)++;
	nnz++;
      }

    // incrementing Ptr
    for (int i = 2; i <= n; i++)
      Ptr(i) += Ptr(i-1);

  }


  //! Conversion from ArrayColSymSparse to symmetric CSC format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& Ind,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();
    int nnz = A.GetDataSize();

    Ptr.Reallocate(n+1);
    Ind.Reallocate(nnz);
    Value.Reallocate(nnz);

    Ptr(0) = 0;
    for (int i = 1; i <= n; i++)
      Ptr(i) = Ptr(i-1) + A.GetColumnSize(i);

    int nb = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
	{
	  Ind(nb) = A.Index(i, j);
	  Value(nb) = A.Value(i, j);
	  nb++;
	}
  }


  //! Conversion from ArrayColSymSparse to CSC format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();

    Vector<Tint, VectFull, Alloc3> IndCol;

    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);

    // sorting by columns
    Sort(IndCol, IndRow, Value);

    Ptr.Reallocate(n+1);
    Ptr.Zero();
    // counting number of non-zero entries
    int nnz = 0;
    for (int i = 0; i < IndCol.GetM(); i++)
      {
	Ptr(IndCol(i) + 1)++;
	nnz++;
      }

    // incrementing Ptr
    for (int i = 2; i <= n; i++)
      Ptr(i) += Ptr(i-1);

  }


  //! Conversion from RowSymSparse to symmetric CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    Vector<Tint, VectFull, Alloc3> IndCol;

    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, false);

    // sorting by columns
    Sort(IndCol, IndRow, Value);

    int n = A.GetN();
    int nnz = A.GetDataSize();

    // creating pointer array
    Ptr.Reallocate(n+1);
    Ptr.Fill(0);
    for (int i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;

    for (int i = 0; i < n; i++)
      Ptr(i+1) += Ptr(i);
  }


  //! Conversion from RowSymSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();

    Vector<Tint, VectFull, Alloc2> IndCol;

    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);

    // sorting by columns
    Sort(IndCol, IndRow, Value);

    Ptr.Reallocate(n+1);
    Ptr.Zero();
    // counting number of non-zero entries
    int nnz = 0;
    for (int i = 0; i < IndCol.GetM(); i++)
      {
	Ptr(IndCol(i) + 1)++;
	nnz++;
      }

    // incrementing Ptr
    for (int i = 2; i <= n; i++)
      Ptr(i) += Ptr(i-1);

  }


  //! Conversion from ArrayRowSymSparse to symmetric CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    Vector<Tint, VectFull, Alloc3> IndCol;

    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, false);

    // sorting by columns
    Sort(IndCol, IndRow, Value);

    int n = A.GetN();
    int nnz = A.GetDataSize();

    // creating pointer array
    Ptr.Reallocate(n+1);
    Ptr.Fill(0);
    for (int i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;

    for (int i = 0; i < n; i++)
      Ptr(i+1) += Ptr(i);
  }


  //! Conversion from ArrayRowSymSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& Ind,
                    Vector<T, VectFull, Alloc4>& AllVal, bool sym_pat)
  {
    int i, j;

    int nnz = A.GetDataSize();
    int n = A.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    Vector<T> Val(nnz);
    int ind = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) != i)
	  {
	    IndRow(ind) = i;
	    IndCol(ind) = A.Index(i, j);
	    Val(ind) = A.Value(i, j);
	    ind++;
	  }

    Sort(ind, IndCol, IndRow, Val);

    Ptr.Reallocate(n+1);
    Ind.Reallocate(nnz + ind);
    AllVal.Reallocate(nnz+ind);
    nnz = ind;
    ind = 0;

    int offset = 0; Ptr(0) = 0;
    for (i = 0; i < n; i++)
      {
	int first_index = ind;
	while (ind < nnz && IndCol(ind) <= i)
	  ind++;

        int size_lower = ind - first_index;
	int size_upper = A.GetRowSize(i);
	int size_row = size_lower + size_upper;

	ind = first_index;
	for (j = 0; j < size_lower; j++)
	  {
	    Ind(offset+j) = IndRow(ind);
	    AllVal(offset+j) = Val(ind);
            ind++;
	  }

	for (j = 0; j < size_upper; j++)
	  {
	    Ind(offset + size_lower + j) = A.Index(i, j);
	    AllVal(offset + size_lower + j) = A.Value(i, j);
          }

        offset += size_row; Ptr(i+1) = offset;
      }
  }


  //! B = A.
  template<class T, class Prop, class Storage, class Allocator>
  void Copy(const Matrix<T, Prop, Storage, Allocator>& A,
	    Matrix<T, Prop, Storage, Allocator>& B)
  {
    B = A;
  }


  //! Conversion from ArrayColSparse to ColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& mat_array,
	    Matrix<T1, Prop1, ColSparse, Allocator1>& mat_csc)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<int, VectFull, CallocAlloc<int> > IndRow;
    Vector<int, VectFull, CallocAlloc<int> > IndCol;

    General sym;
    ConvertToCSC(mat_array, sym, IndCol, IndRow, Val);

    int m = mat_array.GetM();
    int n = mat_array.GetN();

    mat_csc.SetData(m, n, Val, IndCol, IndRow);
  }


  //! Conversion from row-sparse to column-sparse.
  template<class T, class Prop, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop, RowSparse, Alloc1>& A,
	    Matrix<T, Prop, ColSparse, Alloc2>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr;
    Vector<int, VectFull, CallocAlloc<int> > Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }


  //! Conversion from ArrayRowSparse to ColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, ColSparse, Allocator1>& mat_csr)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Val;

    int m = mat_array.GetM();
    int n = mat_array.GetN();
    General sym;
    ConvertToCSC(mat_array, sym, Ptr, IndRow, Val);

    mat_csr.SetData(m, n, Val, Ptr, IndRow);
  }


  //! Conversion from RowSymSparse to ColSparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, RowSymSparse, Alloc1>& A,
	    Matrix<T, Prop2, ColSparse, Alloc2>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr;
    Vector<int, VectFull, CallocAlloc<int> > Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr, Ind;
    Vector<T1> AllVal;

    int n = A.GetM();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, AllVal);

    B.SetData(n, n, AllVal, Ptr, Ind);
  }


  //! Conversion from ColSymSparse to ColSparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, ColSymSparse, Alloc1>& A,
	    Matrix<T, Prop2, ColSparse, Alloc2>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr;
    Vector<int, VectFull, CallocAlloc<int> > Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }


  //! Conversion from ArrayColSymSparse to ColSparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, ArrayColSymSparse, Alloc1>& A,
	    Matrix<T, Prop2, ColSparse, Alloc2>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr;
    Vector<int, VectFull, CallocAlloc<int> > Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }


  /*
    From Sparse formats to CSR format
  */


  //! Conversion from ColSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    int m = A.GetM();
    int n = A.GetN();
    int  nnz = A.GetDataSize();
    if (m <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }

    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();

    // Computation of the indexes of the beginning of rows.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    // Counting the number of entries per row.
    for (int i = 0; i < nnz; i++)
      Ptr(ind_[i])++;

    // Incrementing in order to get the row indexes.
    int increment = 0, size, num_row;
    for (int i = 0; i < m; i++)
      {
	size = Ptr(i);
	Ptr(i) = increment;
	increment += size;
      }

    // Last index.
    Ptr(m) = increment;

    // 'Offset' will be used to get current positions of new entries.
    Vector<int, VectFull, CallocAlloc<int> > Offset(Ptr);
    IndCol.Reallocate(nnz);
    Value.Reallocate(nnz);

    // Loop over the columns.
    for (int j = 0; j < n; j++)
      for (int i = ptr_[j]; i < ptr_[j + 1]; i++)
	{
	  num_row = ind_[i];
	  IndCol(Offset(num_row)) = j;
	  Value(Offset(num_row)) = data_[i];
	  Offset(num_row)++;
	}
  }


  //! Conversion from ArrayColSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    int m = A.GetM();
    int n = A.GetN();
    int  nnz = A.GetDataSize();
    if (m <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }

    // Computation of the indexes of the beginning of rows.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    // Counting the number of entries per row.
    for (int i = 0; i < n; i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
	Ptr(A.Index(i, j))++;

    // Incrementing in order to get the row indexes.
    int increment = 0, size, num_row;
    for (int i = 0; i < m; i++)
      {
	size = Ptr(i);
	Ptr(i) = increment;
	increment += size;
      }

    // Last index.
    Ptr(m) = increment;

    // 'Offset' will be used to get current positions of new entries.
    Vector<int, VectFull, CallocAlloc<int> > Offset(Ptr);
    IndCol.Reallocate(nnz);
    Value.Reallocate(nnz);

    // Loop over the columns.
    for (int j = 0; j < n; j++)
      for (int i = 0; i < A.GetColumnSize(j); i++)
	{
	  num_row = A.Index(j, i);
	  IndCol(Offset(num_row)) = j;
	  Value(Offset(num_row)) = A.Value(j, i);
	  Offset(num_row)++;
	}
  }


  //! Conversion from ColSymSparse to symmetric CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    Vector<Tint, VectFull, Alloc3> IndRow;

    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, false);

    // sorting by rows
    Sort(IndRow, IndCol, Value);

    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;

    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
  }


  //! Conversion from ColSymSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    Vector<Tint, VectFull, Alloc3> IndRow;

    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);

    // sorting by rows
    Sort(IndRow, IndCol, Value);

    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;

    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);

  }


  //! Conversion from ArrayColSymSparse to symmetric CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    Vector<Tint, VectFull, Alloc3> IndRow;

    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, false);

    // sorting by rows
    Sort(IndRow, IndCol, Value);

    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;

    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
  }


  //! Conversion from ArrayColSymSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    Vector<Tint, VectFull, Alloc3> IndRow;

    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);

    // sorting by rows
    Sort(IndRow, IndCol, Value);

    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;

    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
  }


  //! Conversion from ArrayRowSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& IndRow,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Matrix (m,n) with 'nnz' entries.
    int nnz = A.GetDataSize();
    int m = A.GetM();

    // Allocating arrays needed for CSR format.
    Val.Reallocate(nnz);
    IndRow.Reallocate(m + 1);
    IndCol.Reallocate(nnz);

    // Filling the arrays.
    int ind = 0;
    IndRow(0) = 0;
    for (int i = 0; i < m; i++)
      {
	for (int k = 0; k < A.GetRowSize(i); k++)
	  {
	    IndCol(ind) = A.Index(i, k);
	    Val(ind) = A.Value(i, k);
	    ind++;
	  }
	IndRow(i + 1) = ind;
      }
  }


  //! Conversion from RowSymSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& IndRow,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Number of rows and non-zero entries.
    int nnz = A.GetDataSize();
    int m = A.GetM();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();

    // Allocation of arrays for CSR format.
    Val.Reallocate(nnz);
    IndRow.Reallocate(m + 1);
    IndCol.Reallocate(nnz);

    int ind = 0;
    IndRow(0) = 0;
    for (int i = 0; i < m; i++)
      {
	for (int k = ptr_[i]; k < ptr_[i+1]; k++)
	  {
	    IndCol(ind) = ind_[k];
	    Val(ind) = data_[k];
	    ind++;
	  }

	IndRow(i + 1) = ind;
      }
  }


  //! Conversion from ArrayRowSymSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& IndRow,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Number of rows and non-zero entries.
    int nnz = A.GetDataSize();
    int m = A.GetM();

    // Allocation of arrays for CSR format.
    Val.Reallocate(nnz);
    IndRow.Reallocate(m + 1);
    IndCol.Reallocate(nnz);

    int ind = 0;
    IndRow(0) = 0;
    for (int i = 0; i < m; i++)
      {
	for (int k = 0; k < A.GetRowSize(i); k++)
	  {
	    IndCol(ind) = A.Index(i, k);
	    Val(ind) = A.Value(i, k);
	    ind++;
	  }
	IndRow(i + 1) = ind;
      }
  }


  //! Conversion from ColSymSparse to RowSymSparse
  template<class T, class Prop, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
	    Matrix<T, Prop, RowSymSparse, Alloc2>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr, Ind;
    Vector<T, VectFull, Alloc2> Value;

    Symmetric sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }


  //! Conversion from column-oriented sparse to row-oriented sparse.
  /*!
    \param[in] A matrix to be converted.
    \param[out] B converted matrix.
  */
  template<class T1, class T2, class Prop1, class Prop2,
           class Alloc1, class Alloc2>
  void Copy(const Matrix<T1, Prop1, ColSparse, Alloc1>& A,
	    Matrix<T2, Prop2, RowSparse, Alloc2>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr, Ind;
    Vector<T1, VectFull, Alloc2> Value;

    General sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }


  //! Conversion from ArrayColSparse to RowSparse
  /*!
    \param[in] A matrix to be converted.
    \param[out] B converted matrix.
  */
  template<class T1, class T2, class Prop1, class Prop2,
           class Alloc1, class Alloc2>
  void Copy(const Matrix<T1, Prop1, ArrayColSparse, Alloc1>& A,
	    Matrix<T2, Prop2, RowSparse, Alloc2>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr, Ind;
    Vector<T1, VectFull, Alloc2> Value;

    General sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }


  //! Conversion from ArrayRowSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<int, VectFull, CallocAlloc<int> > IndRow;
    Vector<int, VectFull, CallocAlloc<int> > IndCol;

    General sym;
    ConvertToCSR(mat_array, sym, IndRow, IndCol, Val);

    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Val, IndRow, IndCol);
  }


  //! Conversion from ArrayRowSymSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<int, VectFull, CallocAlloc<int> > IndRow;
    Vector<int, VectFull, CallocAlloc<int> > IndCol;

    Symmetric sym;
    ConvertToCSR(mat_array, sym, IndRow, IndCol, Val);

    int m = mat_array.GetM();
    mat_csr.SetData(m, m, Val, IndRow, IndCol);
  }


  /******************************
   * From Sparse to ArraySparse *
   ******************************/


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    int n = A.GetM();
    if (n <= 0)
      {
	B.Clear();
	return;
      }

    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();

    B.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
	int size_row = ptr_[i+1] - ptr_[i];
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = ind_[ptr_[i] + j];
	    B.Value(i, j) = data_[ptr_[i] + j];
	  }
      }

  }


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ColSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    int n = A.GetM();
    if (n <= 0)
      {
	B.Clear();
	return;
      }

    Vector<int, VectFull, CallocAlloc<int> > Ptr, Ind;
    Vector<T0, VectFull, Allocator0> Value;

    Symmetric sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);

    B.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
	int size_row = Ptr(i+1) - Ptr(i);
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = Ind(Ptr(i) + j);
	    B.Value(i, j) = Value(Ptr(i) + j);
	  }
      }
  }


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    if (n <= 0)
      {
	B.Clear();
	return;
      }

    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();

    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
	int size_row = ptr_[i+1] - ptr_[i];
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = ind_[ptr_[i] + j];
	    B.Value(i, j) = data_[ptr_[i] + j];
	  }
      }

  }


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ColSparse, Allocator0>& Acsc,
	    Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int m = Acsc.GetM();
    int n = Acsc.GetN();
    if (n <= 0)
      {
	B.Clear();
	return;
      }

    // conversion to RowSparse
    Matrix<T0, Prop0, RowSparse, Allocator0> A;
    Copy(Acsc, A);

    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();

    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
	int size_row = ptr_[i+1] - ptr_[i];
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = ind_[ptr_[i] + j];
	    B.Value(i, j) = data_[ptr_[i] + j];
	  }
      }
  }


  /***********************************
   * From ArraySparse to ArraySparse *
   ***********************************/


  //! From ArrayRowSymSparse to ArrayRowSymSparse (T0 and T1 different)
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
        int size_row = A.GetRowSize(i);
        B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = A.Index(i, j);
            B.Value(i, j) = A.Value(i, j);
          }
      }
  }


  template<class T, class Prop, class Allocator>
  void Copy(const Matrix<T, Prop, ArrayRowSparse, Allocator>& A,
	    Matrix<T, Prop, ArrayRowSparse, Allocator>& B)
  {
    B = A;
  }


  template<class T, class Prop, class Allocator>
  void Copy(const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A,
	    Matrix<T, Prop, ArrayRowSymSparse, Allocator>& B)
  {
    B = A;
  }


  //! From ArrayRowSparse to ArrayRowSparse (T0 and T1 different)
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
        int size_row = A.GetRowSize(i);
        B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = A.Index(i, j);
            B.Value(i, j) = A.Value(i, j);
          }
      }
  }


  //! do not use
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    abort();
  }


  //! conversion from ArrayColSymSparse to ArrayRowSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayColSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    int n = A.GetM();
    if (n <= 0)
      {
	B.Clear();
	return;
      }

    Vector<int, VectFull, CallocAlloc<int> > Ptr, Ind;
    Vector<T0, VectFull, Allocator0> Value;

    Symmetric sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);

    B.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
	int size_row = Ptr(i+1) - Ptr(i);
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = Ind(Ptr(i) + j);
	    B.Value(i, j) = Value(Ptr(i) + j);
	  }
      }
  }


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& Acsc,
	    Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int m = Acsc.GetM();
    int n = Acsc.GetN();
    if (n <= 0)
      {
	B.Clear();
	return;
      }

    // conversion to RowSparse
    Matrix<T0, Prop0, RowSparse, Allocator0> A;
    Copy(Acsc, A);

    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();

    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
	int size_row = ptr_[i+1] - ptr_[i];
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = ind_[ptr_[i] + j];
	    B.Value(i, j) = data_[ptr_[i] + j];
	  }
      }
  }


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int i, j;

    int nnz = A.GetDataSize();
    int n = A.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz),IndCol(nnz);
    Vector<T1, VectFull, Allocator1> Val(nnz);
    int ind = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) != i)
	  {
	    IndRow(ind) = i;
	    IndCol(ind) = A.Index(i, j);
	    Val(ind) = A.Value(i, j);
	    ind++;
	  }
    Sort(ind, IndCol, IndRow, Val);
    nnz = ind;
    ind = 0;

    B.Reallocate(n, n);
    for (i = 0; i < n; i++)
      {
	int first_index = ind;
	while (ind < nnz && IndCol(ind) <= i)
	  ind++;
	int size_lower = ind - first_index;
	int size_upper = A.GetRowSize(i);
	int size_row = size_lower + size_upper;
	B.ResizeRow(i, size_row);
	ind = first_index;
	for (j = 0; j < size_lower; j++)
	  {
	    B.Index(i, j) = IndRow(ind);
	    B.Value(i, j) = Val(ind);
	    ind++;
	  }
	for (j = 0; j < size_upper; j++)
	  {
	    B.Index(i, size_lower + j) = A.Index(i, j);
	    B.Value(i, size_lower + j) = A.Value(i, j);
	  }
	B.AssembleRow(i);
      }
  }


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B)
  {
    int i, j;

    int nnz = A.GetDataSize();
    int n = A.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    Vector<T1> Val(nnz);
    int ind = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) != i)
	  {
	    IndRow(ind) = i;
	    IndCol(ind) = A.Index(i, j);
	    Val(ind) = A.Value(i, j);
	    ind++;
	  }
    Sort(ind, IndCol, IndRow, Val);
    nnz = ind;
    ind = 0;

    B.Reallocate(n, n);
    for (i = 0; i < n; i++)
      {
	int first_index = ind;
	while (ind < nnz && IndCol(ind) <= i)
	  ind++;
	int size_lower = ind - first_index;
	int size_upper = A.GetRowSize(i);
	int size_row = size_lower + size_upper;
	B.ResizeColumn(i, size_row);
	ind = first_index;
	for (j = 0; j < size_lower; j++)
	  {
	    B.Index(i, j) = IndRow(ind);
	    B.Value(i, j) = Val(ind);
	    ind++;
	  }
	for (j = 0; j < size_upper; j++)
	  {
	    B.Index(i, size_lower + j) = A.Index(i, j);
	    B.Value(i, size_lower + j) = A.Value(i, j);
	  }
	B.AssembleColumn(i);
      }
  }


  //! Conversion from ArrayRowSparse to ArrayColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& A,
	    Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B)
  {
    int i;

    // Matrix (m,n) with nnz entries.
    int nnz = A.GetDataSize();
    int n = A.GetN();

    // Conversion in coordinate format.
    Vector<T1> Val;
    Vector<int, VectFull, CallocAlloc<int> > IndRow, IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);

    // Sorting with respect to column numbers.
    Sort(IndCol, IndRow, Val);

    // Constructing pointer array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(n + 1);

    // Counting non-zero entries per column.
    for (i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);

    // we fill matrix B
    for (int i = 0; i < n; i++)
      {
	int size_col = Ptr(i+1) - Ptr(i);
	if (size_col > 0)
	  {
	    B.ReallocateColumn(i, size_col);
	    for (int j = Ptr(i); j < Ptr(i+1); j++)
	      {
		B.Index(i, j-Ptr(i)) = IndRow(j);
		B.Value(i, j-Ptr(i)) = Val(j);
	      }
	  }
      }
  }


  /**************************
   * ComplexSparse matrices *
   **************************/


  //! Conversion from ArrayRowComplexSparse to RowComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  Copy(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowComplexSparse, Allocator1>& mat_csr)
  {
    int i, k;

    // Non-zero entries (real and imaginary part).
    int nnz_real = mat_array.GetRealDataSize();
    int nnz_imag = mat_array.GetImagDataSize();
    int m = mat_array.GetM();

    // Allocation of arrays for CSR.
    Vector<T0, VectFull, Allocator1> Val_real(nnz_real), Val_imag(nnz_imag);
    Vector<int, VectFull, CallocAlloc<int> > IndRow_real(m + 1);
    Vector<int, VectFull, CallocAlloc<int> > IndRow_imag(m + 1);
    Vector<int, VectFull, CallocAlloc<int> > IndCol_real(nnz_real);
    Vector<int, VectFull, CallocAlloc<int> > IndCol_imag(nnz_imag);

    int ind_real = 0, ind_imag = 0;
    IndRow_real(0) = 0;
    IndRow_imag(0) = 0;
    // Loop over rows.
    for (i = 0; i < m; i++)
      {
	for (k = 0; k < mat_array.GetRealRowSize(i); k++)
	  {
	    IndCol_real(ind_real) = mat_array.IndexReal(i, k);
	    Val_real(ind_real) = mat_array.ValueReal(i, k);
	    ind_real++;
	  }

	IndRow_real(i + 1) = ind_real;
	for (k = 0; k < mat_array.GetImagRowSize(i); k++)
	  {
	    IndCol_imag(ind_imag) = mat_array.IndexImag(i, k);
	    Val_imag(ind_imag) = mat_array.ValueImag(i, k);
	    ind_imag++;
	  }

	IndRow_imag(i + 1) = ind_imag;
      }

    mat_csr.SetData(m, m, Val_real, IndRow_real, IndCol_real,
		    Val_imag, IndRow_imag, IndCol_imag);
  }


  //! Conversion from ArrayRowSymComplexSparse to RowSymComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  Copy(const Matrix<T0, Prop0,
       ArrayRowSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& mat_csr)
  {
    int i, k;

    // Non-zero entries (real and imaginary part).
    int nnz_real = mat_array.GetRealDataSize();
    int nnz_imag = mat_array.GetImagDataSize();
    int m = mat_array.GetM();

    // Allocation of arrays for CSR.
    Vector<T0, VectFull, Allocator1> Val_real(nnz_real), Val_imag(nnz_imag);
    Vector<int, VectFull, CallocAlloc<int> > IndRow_real(m + 1);
    Vector<int, VectFull, CallocAlloc<int> > IndRow_imag(m + 1);
    Vector<int, VectFull, CallocAlloc<int> > IndCol_real(nnz_real);
    Vector<int, VectFull, CallocAlloc<int> > IndCol_imag(nnz_imag);

    int ind_real = 0, ind_imag = 0;
    IndRow_real(0) = 0;
    IndRow_imag(0) = 0;
    // Loop over rows.
    for (i = 0; i < m; i++)
      {
	for (k = 0; k < mat_array.GetRealRowSize(i); k++)
	  {
	    IndCol_real(ind_real) = mat_array.IndexReal(i, k);
	    Val_real(ind_real) = mat_array.ValueReal(i, k);
	    ind_real++;
	  }

	IndRow_real(i + 1) = ind_real;
	for (int k = 0; k < mat_array.GetImagRowSize(i); k++)
	  {
	    IndCol_imag(ind_imag) = mat_array.IndexImag(i, k);
	    Val_imag(ind_imag) = mat_array.ValueImag(i, k);
	    ind_imag++;
	  }

	IndRow_imag(i + 1) = ind_imag;
      }

    mat_csr.SetData(m, m, Val_real, IndRow_real, IndCol_real,
		    Val_imag, IndRow_imag, IndCol_imag);
  }


  /***********************
   * GetSymmetricPattern *
   ***********************/


  template<class T, class Prop, class Storage, class Allocator,
           class Tint, class Allocator2, class Allocator3>
  void GetSymmetricPattern(const Matrix<T, Prop, Storage, Allocator>& A,
                           Vector<Tint, VectFull, Allocator2>& Ptr,
                           Vector<Tint, VectFull, Allocator3>& Ind)
  {
    int n = A.GetM();

    // Converting to coordinates.
    Vector<int, VectFull, CallocAlloc<int> > IndRow, IndCol;
    Vector<T> Value;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value);

    // clearing values
    Value.Clear();

    // Sorting columns too.
    Vector<int, VectFull, CallocAlloc<int> > IndRow2, IndCol2, Index(2*n);
    IndRow2 = IndRow;
    IndCol2 = IndCol;
    Sort(IndCol2.GetM(), IndCol2, IndRow2);

    Tint max_nnz = 0;
    for (int i = 0; i < IndRow.GetM(); i++)
      if (IndRow(i) <= IndCol(i))
        max_nnz++;

    for (int i = 0; i < IndRow.GetM(); i++)
      if (IndCol2(i) <= IndRow2(i))
        max_nnz++;

    // then symmetrization of pattern and conversion to csr.
    Ptr.Reallocate(n+1);
    Ind.Reallocate(max_nnz);
    Tint j_begin = 0, j_end = 0;
    int size_row = 0;
    Tint j2_begin = 0, j2_end = 0;
    Ptr(0) = 0;
    for (int i = 0; i < A.GetM(); i++)
      {
        j_begin = j_end;
        size_row = 0;
        // We retrieve column numbers.
        while ( (j_end < IndRow.GetM()) && (IndRow(j_end) == i))
          {
            if (IndRow(j_end) <= IndCol(j_end))
              {
                Index(size_row) = IndCol(j_end);
                size_row++;
              }

            j_end++;
          }

        j2_begin = j2_end;
        while ( (j2_end < IndCol2.GetM()) && (IndCol2(j2_end) == i))
          {
            if (IndCol2(j2_end) <= IndRow2(j2_end))
              {
                Index(size_row) = IndRow2(j2_end);
                size_row++;
              }

            j2_end++;
          }

        // Sorting indexes.
        Assemble(size_row, Index);

        // Updating Ptr, Ind.
        for (int j = 0; j < size_row; j++)
	  Ind(Ptr(i) + j) = Index(j);

        Ptr(i+1) = Ptr(i) + size_row;
      }

    IndRow2.Clear(); IndCol2.Clear();
    IndRow.Clear(); IndCol.Clear();
    Ind.Resize(Ptr(n));
  }


  template<class T, class Prop, class Storage, class Allocator, class AllocI>
  void GetSymmetricPattern(const Matrix<T, Prop, Storage, Allocator>& A,
                           Matrix<int, Symmetric, RowSymSparse, AllocI>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr, Ind;

    GetSymmetricPattern(A, Ptr, Ind);

    int n = A.GetM();
    Vector<int, VectFull, CallocAlloc<int> > Val(Ptr(n));
    // We put Ptr and Ind into the matrix B.
    B.SetData(n, n, Val, Ptr, Ind);
  }


  /*****************************************************
   * Conversion from sparse matrices to dense matrices *
   *****************************************************/


  //! Conversion from RowSparse to RowMajor.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void Copy(Matrix<T, Prop, RowSparse, Allocator1>& A,
            Matrix<T, Prop, RowMajor, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* data = A.GetData();

    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = ptr[i]; j < ptr[i+1]; j++)
        B(i, ind[j]) = data[j];

  }


  //! Conversion from ArrayRowSparse to RowMajor.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void Copy(Matrix<T, Prop, ArrayRowSparse, Allocator1>& A,
            Matrix<T, Prop, RowMajor, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();

    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        B(i, A.Index(i, j)) = A.Value(i, j);

  }


  //! Conversion from ArrayRowSymSparse to RowSymPacked.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void Copy(Matrix<T, Prop, ArrayRowSymSparse, Allocator1>& A,
            Matrix<T, Prop, RowSymPacked, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	B(i, A.Index(i, j)) = A.Value(i, j);

  }


  /*****************************************************
   * Conversion from dense matrices to sparse matrices *
   *****************************************************/


  //! Conversion from RowSymPacked to RowSymSparse.
  template<class T>
  void ConvertToSparse(const Matrix<T, Symmetric, RowSymPacked>& A,
                       Matrix<T, Symmetric, RowSymSparse>& B,
		       const T& threshold)
  {
    int nnz = 0;
    int n = A.GetM();
    for (int i = 0; i < n; i++)
      for (int j = i; j < n; j++)
        if (abs(A(i, j)) > threshold)
          nnz++;

    IVect IndCol(nnz), IndRow(n+1);
    Vector<T> Value(nnz);
    nnz = 0; IndRow(0) = 0;
    for (int i = 0; i < n; i++)
      {
        for (int j = i; j < n; j++)
          if (abs(A(i, j)) > threshold)
            {
              IndCol(nnz) = j;
              Value(nnz) = A(i, j);
              nnz++;
            }

        IndRow(i+1) = nnz;
      }

    B.SetData(n, n, Value, IndRow, IndCol);

  }


  //! Conversion from RowMajor to ArrayRowSparse.
  template<class T>
  void ConvertToSparse(const Matrix<T, General, RowMajor>& A,
                       Matrix<T, General, ArrayRowSparse>& B,
		       const T& threshold)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
        int size_row = 0;
        for (int j = 0; j < n; j++)
          if (abs(A(i, j)) > threshold)
            size_row++;

        B.ReallocateRow(i, size_row);

        size_row = 0;
        for (int j = 0; j < n; j++)
          if (abs(A(i, j)) > threshold)
            {
              B.Index(i, size_row) = j;
              B.Value(i, size_row) = A(i, j);
              size_row++;
            }
      }
  }


  //! Conversion from RowMajor to RowSparse.
  template<class T>
  void ConvertToSparse(const Matrix<T, General, RowMajor>& A,
                       Matrix<T, General, RowSparse>& B,
		       const T& threshold)
  {
    int nnz = 0;
    int m = A.GetM();
    int n = A.GetN();
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        if (abs(A(i, j)) > threshold)
          nnz++;

    Vector<int, VectFull, CallocAlloc<int> > IndCol(nnz), IndRow(m+1);
    Vector<T> Value(nnz);
    nnz = 0; IndRow(0) = 0;
    for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
          if (abs(A(i, j)) > threshold)
            {
              IndCol(nnz) = j;
              Value(nnz) = A(i, j);
              nnz++;
            }

        IndRow(i+1) = nnz;
      }

    B.SetData(m, n, Value, IndRow, IndCol);

  }

} // namespace Seldon.

#define SELDON_FILE_MATRIX_CONVERSIONS_CXX
#endif
