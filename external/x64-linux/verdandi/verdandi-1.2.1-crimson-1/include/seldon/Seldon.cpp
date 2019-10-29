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


#ifndef SELDON_FILE_SELDON_CPP

#include "Seldon.hxx"

namespace Seldon
{
  template class MallocAlloc<int>;
  template class Vector_Base<int, MallocAlloc<int> >;
  template class Vector<int, VectFull, MallocAlloc<int> >;
  template class Vector<Vector<int, VectFull, MallocAlloc<int> >, Collection,
                MallocAlloc<Vector<int, VectFull, MallocAlloc<int> > > >;
  template class MallocAlloc<double>;
  template class Vector_Base<double, MallocAlloc<double> >;
  template class Vector<double, VectFull, MallocAlloc<double> >;
  template class Vector<Vector<double, VectFull, MallocAlloc<double> >,
            Collection, MallocAlloc<Vector<double, VectFull, MallocAlloc<double> > > >;

  template class Matrix_Base<int, MallocAlloc<int> >;
  template class Matrix_Pointers<int, General, RowMajor, MallocAlloc<int> >;
  template class Matrix<int, General, RowMajor, MallocAlloc<int> >;
  template class Matrix_Base<double, MallocAlloc<double> >;
  template class Matrix_Pointers<double, General, RowMajor, MallocAlloc<double> >;
  template class Matrix<double, General, RowMajor, MallocAlloc<double> >;

  template class Vector<double, VectSparse, MallocAlloc<double> >;

  template class Matrix_Sparse<double, General, RowSparse, MallocAlloc<double> >;
  template class Matrix<double, General, RowSparse, MallocAlloc<double> >;

  template void ConvertMatrix_to_Coordinates(const Matrix<double, General, RowSparse, MallocAlloc<double> >& A,
                                             Vector<int>& IndRow, Vector<int>& IndCol,
                                             Vector<double, VectFull>& Val,
                                             int index = 0, bool sym = false);

  template void ConvertMatrix_from_Coordinates(Vector<int>& IndRow, Vector<int>& IndCol,
                                               Vector<double, VectFull>& Val,
                                               Matrix<double, General, RowSparse, MallocAlloc<double> >& A,
                                               int index = 0);

  // Skips one vector in an input stream.
  void skip_vector_double(istream& input_stream)
  {
    int size;
    input_stream.read(reinterpret_cast<char*>(&size), sizeof(int));
    input_stream.seekg(size * sizeof(double), istream::cur);
  }

  // Skips one matrix in an input stream.
  void skip_matrix_double(istream& input_stream)
  {
    int m, n;
    input_stream.read(reinterpret_cast<char*>(&m), sizeof(int));
    input_stream.read(reinterpret_cast<char*>(&n), sizeof(int));
    input_stream.seekg(m * n * sizeof(double), istream::cur);
  }
}


#define SELDON_FILE_SELDON_CPP
#endif
