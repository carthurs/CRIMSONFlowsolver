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


#ifndef SELDON_FILE_MATRIX_SPARSE_IOMATRIXMARKET_HXX


namespace Seldon
{


  template <class Prop, class Allocator>
  void ReadHarwellBoeing(string filename,
                         Matrix<float, Prop, ColSparse, Allocator>& A);


  template <class Prop, class Allocator>
  void ReadHarwellBoeing(string filename,
                         Matrix<double, Prop, ColSparse, Allocator>& A);


  template <class T, class Prop, class Storage, class Allocator>
  void ReadHarwellBoeing(string filename,
                         string value_type, string matrix_type,
                         Matrix<T, Prop, Storage, Allocator> & A);


  template<class T>
  void PrintComplexValuesHarwell(int nnz, const Vector<complex<T> >& Val,
                                 ofstream& file_out);


  template<class T>
  void PrintComplexValuesHarwell(int nnz, const Vector<T>& Val,
                                 ofstream& file_out);


  template<class T, class Prop, class Storage, class Allocator>
  void WriteHarwellBoeing(const Matrix<T, Prop, Storage, Allocator>& A,
                          const string& file_name, bool complex);


  template<class T, class Prop, class Storage, class Allocator>
  void WriteHarwellBoeing(const Matrix<T, Prop, Storage, Allocator>& A,
                          const string& file_name);


  template<class T, class Prop, class Storage, class Allocator>
  void WriteHarwellBoeing(const Matrix<complex<T>,
                          Prop, Storage, Allocator>& A,
                          const string& file_name);

}


#define SELDON_FILE_MATRIX_SPARSE_IOMATRIXMARKET_HXX
#endif
