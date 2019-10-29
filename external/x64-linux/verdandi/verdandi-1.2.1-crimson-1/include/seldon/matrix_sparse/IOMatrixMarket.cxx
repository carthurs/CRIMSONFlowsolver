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


#ifndef SELDON_FILE_MATRIX_SPARSE_IOMATRIXMARKET_CXX


#include "IOMatrixMarket.hxx"
#include <iomanip>

namespace Seldon
{


  //! Reads a sparse matrix in a file in Harwell-Boeing format.
  /*! This functions was written to read the files in Harwell-Boeing format as
    distributed on the Matrix Market, http://math.nist.gov/MatrixMarket/. A
    file from Matrix Market is associated with a type encoded in three
    characters (which is usually its extension). The supported types have:
    - 'R' (real) or 'P' (pattern) as first character; complex-valued matrices
    ('C') are not supported;
    - 'U' (unsymmetric) or 'R' (rectangular) as second character; 'S'
    (symmetric), 'H' (Hermitian) and 'Z' (skew symmetric) are not supported;
    - 'A' (assembled) as third character; 'E' (elemental) is not supported.
    \param[in] filename path to the file that contains the matrix.
    \param[out] A the matrix to be structured and filled with the values of \a
    filename.
  */
  template <class Prop, class Allocator>
  void ReadHarwellBoeing(string filename,
                         Matrix<float, Prop, ColSparse, Allocator>& A)
  {
    ReadHarwellBoeing(filename, "real", "general", A);
  }


  //! \copydoc ReadHarwellBoeing(string filename, Matrix<float, Prop, ColSparse, Allocator>& A)
  template <class Prop, class Allocator>
  void ReadHarwellBoeing(string filename,
                         Matrix<double, Prop, ColSparse, Allocator>& A)
  {
    ReadHarwellBoeing(filename, "real", "general", A);
  }


  //! Reads a sparse matrix in a file in Harwell-Boeing format.
  /*! This functions was written to read the files in Harwell-Boeing format as
    distributed on the Matrix Market, http://math.nist.gov/MatrixMarket/. A
    file from Matrix Market is associated with a type encoded in three
    characters (which is usually its extension). The supported types have:
    - 'R' (real) or 'P' (pattern) as first character; complex-valued matrices
    ('C') are not supported;
    - 'U' (unsymmetric) or 'R' (rectangular) as second character; 'S'
    (symmetric), 'H' (Hermitian) and 'Z' (skew symmetric) are not supported;
    - 'A' (assembled) as third character; 'E' (elemental) is not supported.
    \param[in] filename path to the file that contains the matrix.
    \param[in] value_type the type of involved data: only "real" is supported.
    \param[in] matrix_type the type of matrix: only "general" is supported.
    \param[out] A the matrix to be structured and filled with the values of \a
    filename.
    \note This function is not supposed to be called directly unless the user
    knows what he/she is doing. See overloaded functions ReadHarwellBoeing
    with only \a filename and \a A as arguments.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void ReadHarwellBoeing(string filename,
                         string value_type, string matrix_type,
                         Matrix<T, Prop, Storage, Allocator> & A)
  {
    int i, j, k;
    string line, element;

    // Dimensions of the output matrix.
    int Nrow, Ncol;
    int Nnonzero;

    ifstream input_stream(filename.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!input_stream.good())
      throw IOError("ReadHarwellBoeing(string filename, "
                    "Matrix& A, string value_type, string matrix_type)",
                    "Unable to read file \"" + filename + "\".");
#endif

    /*** Parsing headers ***/

    // First line.
    string title, name;
    getline(input_stream, line);
    title = line.substr(0, 72);
    name = line.substr(72, 8);

    // Second line.
    int Nline_ptr, Nline_ind, Nline_val, Nline_rhs;
    input_stream >> i >> Nline_ptr >> Nline_ind >> Nline_val >> Nline_rhs;

    // Third line.
    string type;
    int Nelemental;
    input_stream >> type >> Nrow >> Ncol >> Nnonzero >> Nelemental;

    if (type.size() != 3)
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains a matrix of "
                          "type \"" + type + "\", which is not a valid type. "
                          "The type must contain exactly three characters.");

    if (type.substr(0, 1) != "R" && type.substr(0, 1) != "C"
        && type.substr(0, 1) != "P")
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains a matrix of "
                          "type \"" + type + "\", which is not a valid type. "
                          "The first character of that type must be 'R', 'C' "
                          "(not supported anyway) or 'P'.");

    if (type.substr(0, 1) == "C")
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains a "
                          "complex-valued matrix, which is not supported.");

    if (type.substr(0, 1) == "R" && value_type != "real")
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains a real-valued "
                          "matrix, while the input matrix 'A' is not "
                          "declared as such.");

    if (type.substr(1, 1) == "H")
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains a Hermitian "
                          "matrix, which is not supported.");

    if (type.substr(1, 1) == "Z")
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains a skew "
                          "symmetric matrix, which is not supported.");

    if (type.substr(1, 1) == "S")
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains a "
                          "symmetric matrix, which is not supported.");

    if (Nelemental != 0 || type.substr(2, 1) == "E")
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains an elemental "
                          "matrix, which is not supported.");

    getline(input_stream, line);

    // Fourth line.
    getline(input_stream, line);
    int line_width_ptr, element_width_ptr, line_width_ind, element_width_ind;
    element = line.substr(0, 16);
    sscanf(element.c_str(), "(%dI%d)", &line_width_ptr, &element_width_ptr);
    element = line.substr(16, 16);
    sscanf(element.c_str(), "(%dI%d)", &line_width_ind, &element_width_ind);
    int line_width_val = 0;
    int element_width_val = 0;
    if (type.substr(0, 1) != "P")
      {
        element = line.substr(32, 16);

        // Splits the format to retrieve the useful widths. This part is
        // tricky because no parsing (in C++) of Fortran format was found.
        vector<int> vect;
        string delimiter = " (ABEFGILOPZ*.)";
        string tmp_str;
        int tmp_int;
        string::size_type index_beg = element.find_first_not_of(delimiter);
        string::size_type index_end;
        while (index_beg != string::npos)
          {
            index_end = element.find_first_of(delimiter, index_beg);
            tmp_str = element.substr(index_beg,
                                     index_end == string::npos ?
                                     string::npos : (index_end-index_beg));
            istringstream(tmp_str) >> tmp_int;
            vect.push_back(tmp_int);
            index_beg = element.find_first_not_of(delimiter, index_end);
          }

        if (vect.size() < 3)
          throw WrongArgument("ReadHarwellBoeing(string filename, Matrix& A, "
                              "string value_type, string matrix_type)",
                              "File \"" + filename + "\" contains values "
                              "written in format \"" + element + "\", which "
                              "could not be parsed.");

        line_width_val = vect[vect.size() - 3];
        element_width_val = vect[vect.size() - 2];
      }

    // Fifth header line, if any: ignored. RHS are not read.
    if (Nline_rhs != 0)
      getline(input_stream, line);

    /*** Allocations ***/

    // Content of output matrix A.
    int* A_ptr;
    int* A_ind;
    T* A_data;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	A_ptr = reinterpret_cast<int*>(calloc(Ncol + 1, sizeof(int)));

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
        A_ptr = NULL;
      }

    if (A_ptr == NULL)
      throw NoMemory("ReadHarwellBoeing(string filename, "
                     "Matrix& A, string value_type, string matrix_type)",
		     "Unable to allocate memory for an array of "
		     + to_str(Ncol + 1) + " integers.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

        // Reallocates 'A_ind' and 'A_data' in order to append the
        // elements of the i-th row of C.
        A_ind = reinterpret_cast<int*>(calloc(Nnonzero, sizeof(int)));
        A_data = reinterpret_cast<T*>
          (A.GetAllocator().allocate(Nnonzero));

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
        A_ind = NULL;
        A_data = NULL;
      }

    if (A_ind == NULL || A_data == NULL)
      throw NoMemory("ReadHarwellBoeing(string filename, "
                     "Matrix& A, string value_type, string matrix_type)",
                     "Unable to allocate memory for an array of "
                     + to_str(Nnonzero) + " integers "
                     "and for an array of "
                     + to_str(sizeof(T) * Nnonzero) + " bytes.");
#endif

    /*** Reads the structure ***/

    int index = 0;
    for (i = 0; i < Nline_ptr; i++)
      {
        getline(input_stream, line);
        k = 0;
        for (j = 0; j < line_width_ptr; j++)
          {
            istringstream(line.substr(k, element_width_ptr)) >> A_ptr[index];
            // The indexes are 1-based, so this corrects it:
            A_ptr[index]--;
            index++;
            if (index == Ncol + 1)
              // So as not to read more elements than actually available on
              // the line.
              break;
            k += element_width_ptr;
          }
        if (index == Ncol + 1)
          break;
      }

    index = 0;
    for (i = 0; i < Nline_ind; i++)
      {
        getline(input_stream, line);
        k = 0;
        for (j = 0; j < line_width_ind; j++)
          {
            istringstream(line.substr(k, element_width_ind)) >> A_ind[index];
            // The indexes are 1-based, so this corrects it:
            A_ind[index]--;
            index++;
            if (index == Nnonzero)
              // So as not to read more elements than actually available on
              // the line.
              break;
            k += element_width_ind;
          }
        if (index == Nnonzero)
          break;
      }

    /*** Reads the values ***/

    if (type.substr(0, 1) == "P")
      for (i = 0; i < Nnonzero; i++)
        A_data[i] = T(1);
    else
      {
        index = 0;
        for (i = 0; i < Nline_val; i++)
          {
            getline(input_stream, line);
            k = 0;
            for (j = 0; j < line_width_val; j++)
              {
                istringstream(line.substr(k, element_width_val))
                  >> A_data[index];
                index++;
                if (index == Nnonzero)
                  // So as not to read more elements than actually available
                  // on the line.
                  break;
                k += element_width_val;
              }
            if (index == Nnonzero)
              break;
          }
      }

#ifdef SELDON_CHECK_IO
    // Checks if the stream is still in good shape.
    if (!input_stream.good())
      throw IOError("ReadHarwellBoeing(string filename, "
                    "Matrix& A, string value_type, string matrix_type)",
                    "Unable to read all values in file \""
                    + filename + "\".");
#endif

    input_stream.close();

    A.SetData(Nrow, Ncol, Nnonzero, A_data, A_ptr, A_ind);
  }


  template<class T>
  void PrintComplexValuesHarwell(int nnz, const Vector<complex<T> >& Val,
                                 ofstream& file_out)
  {
    for (int i = 0; i < 2*nnz; i += 3)
      {
        for (int j = i; j < min(i+3, 2*nnz); j++)
          {
            if (j%2 == 0)
              file_out << setw(23) << real(Val(j/2));
            else
              file_out << setw(23) << imag(Val(j/2));
          }

        file_out << '\n';
      }
  }


  template<class T>
  void PrintComplexValuesHarwell(int nnz, const Vector<T>& Val,
                                 ofstream& file_out)
  {
    for (int i = 0; i < nnz; i += 3)
      {
        for (int j = i; j < min(i+3, nnz); j++)
          file_out << setw(23) << Val(j);

        file_out << '\n';
      }
  }


  template<class T, class Prop, class Storage, class Allocator>
  void WriteHarwellBoeing(const Matrix<T, Prop, Storage, Allocator>& A,
                          const string& file_name, bool complex)
  {
    // converting to CSR
    IVect Ptr, Ind;
    Vector<T> Val;
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    // number of columns and non-zero entries
    int nnz = Val.GetM();
    int m = A.GetM();
    int n = A.GetN();

    // First line
    ofstream file_out(file_name.data());

    string title("Seldon");
    string key("S0000008");
    file_out.setf(ios::left);
    file_out << setw(72) << title;
    file_out << setw(8) << key;
    file_out << endl;

    // Compute column pointer format
    int ptr_len = int(ceil( log10(0.1 + nnz + 1) )) + 1;
    int ptr_nperline = min(80 / ptr_len, n+1);
    int ptrcrd = n / ptr_nperline + 1;
    string ptrfmt = string("(") + to_str(ptr_nperline)
      + "I" + to_str(ptr_len) + ")";

    // Compute row index format
    int ind_len = int(ceil( log10(0.1 + n + 1) )) + 1;
    int ind_nperline = min(80 / ind_len, nnz);
    int indcrd = (nnz-1) / ind_nperline + 1;
    string indfmt = string("(") + to_str(ind_nperline)
      + 'I' + to_str(ind_len) + ')';

    // compute value format
    string valfmt("(3D23.15)");
    int valcrd = (nnz-1) / 3 + 1;
    if (complex)
      valcrd = (2*nnz-1) / 3 + 1;

    int rhscrd = 0;
    int totcrd = ptrcrd + indcrd + valcrd + rhscrd;

    // Second line
    file_out << setw(14) << totcrd;
    file_out << setw(14) << ptrcrd;
    file_out << setw(14) << indcrd;
    file_out << setw(14) << valcrd;
    file_out << setw(14) << rhscrd;
    file_out << endl;

    // Third line
    int neltvl = 0;
    char mxtype[4];
    if (complex)
      mxtype[0] = 'C';
    else
      mxtype[0] = 'R';

    if (m == n)
      mxtype[1] = 'U';
    else
      mxtype[1] = 'R';

    mxtype[2] = 'A';
    mxtype[3] = '\0';

    file_out << mxtype << "           ";
    file_out << setw(14) << m;
    file_out << setw(14) << n;
    file_out << setw(14) << nnz;
    file_out << setw(14) << neltvl;
    file_out << endl;

    // Fourth line
    file_out << setw(16) << ptrfmt;
    file_out << setw(16) << indfmt;
    file_out << setw(20) << valfmt;
    file_out << setw(20) << valfmt;
    file_out << endl;

    // printing ptr values
    file_out.setf(ios::left);
    for (int i = 0; i <= n; i += ptr_nperline)
      {
        for (int j = i; j < min(i+ptr_nperline, n+1); j++)
          file_out << setw(ptr_len) << Ptr(j)+1;

        file_out << '\n';
      }

    // printing ind values
    for (int i = 0; i < nnz; i += ind_nperline)
      {
        for (int j = i; j < min(i+ind_nperline, nnz); j++)
          file_out << setw(ind_len) << Ind(j)+1;

        file_out << '\n';
      }

    // printing values
    file_out.precision(15);
    file_out.setf(ios::scientific);
    PrintComplexValuesHarwell(nnz, Val, file_out);

    file_out.close();
  }


  //! A is written in Harwell-Boeing file (.rua extension)
  template<class T, class Prop, class Storage, class Allocator>
  void WriteHarwellBoeing(const Matrix<T, Prop, Storage, Allocator>& A,
                          const string& file_name)
  {
    WriteHarwellBoeing(A, file_name, false);
  }


  //! A is written in Harwell-Boeing file (.cua extension)
  template<class T, class Prop, class Storage, class Allocator>
  void WriteHarwellBoeing(const Matrix<complex<T>,
                          Prop, Storage, Allocator>& A,
                          const string& file_name)
  {
    WriteHarwellBoeing(A, file_name, true);
  }

}


#define SELDON_FILE_MATRIX_SPARSE_IOMATRIXMARKET_CXX
#endif
