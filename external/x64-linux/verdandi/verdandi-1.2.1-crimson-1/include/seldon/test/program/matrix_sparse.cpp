// Copyright (C) 2001-2009 Vivien Mallet
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


#define SELDON_WITH_DEBUG_LEVEL_4

#include "Seldon.hxx"
#include "SeldonSolver.hxx"

using namespace Seldon;

int main(int argc, char** argv)
{

  {
    // We are testing unsymmetric matrix.
    int m = 5, n = 4;
    Matrix<double, General, ArrayRowSparse> A(m,n);

    cout << "Number of rows in the matrix " << A.GetM() << endl;
    cout << "Number of columns in the matrix " << A.GetN() << endl;

    // We can add elements by using the method index and value:
    int irow = 3, icol, nb_values_row = 2;
    A.ReallocateRow(irow, nb_values_row);
    A.Index(irow, 0) = 3; A.Value(irow, 0) = -2.0;
    A.Index(irow, 1) = 1; A.Value(irow, 1) = 3.0;
    // values not sorted you have to call AssembleRow
    A.AssembleRow(irow); DISP(A);

    // with reallocate, A is empty with m columns, n rows
    A.Reallocate(m,n); DISP(A);

    irow = 2;
    A.ReallocateRow(irow, nb_values_row);
    A.Index(irow, 0) = 0; A.Value(irow, 0) = 1.0;
    A.Index(irow, 1) = 2; A.Value(irow, 1) = 1.0;
    // Values already sorted, no need to call AssembleRow.
    DISP(A);


    // We can add elements by using the method AddInteraction.
    double value = 2.5; irow = 1; icol = 3;
    A.AddInteraction(irow, icol, value);
    DISP(A);

    // AddInteractionsRow (no need to sort numbers).
    IVect columns(nb_values_row); Vector<double> values(nb_values_row);
    columns(0) = 3; columns(1) = 2;
    values(0) = -3.0; values(1) = 4.6;
    A.AddInteractionRow(irow, nb_values_row, columns, values);
    DISP(A);

    // AddInteractionsColumn (no need to sort numbers).
    int nb_values_col = 3; icol = 0;
    IVect rows(nb_values_col); values.Reallocate(nb_values_col);
    rows(0) = 3; rows(1) = 0; rows(2) = 1;
    values(0) = -1.0; values(1) = 1.5; values(2) = 0.5;
    A.AddInteractionColumn(icol, nb_values_col, rows, values);
    DISP(A);

    // You can directly access A(i,j) and modify it.
    DISP(A(4,2));
    A.Get(4,2) = 2.5; DISP(A(4,2));

    // You can clear a row,
    A.ClearRow(3);
    // and resize a row and add new coefficients at the end.
    int size_row = A.GetRowSize(0);
    A.ResizeRow(0, size_row+1);
    A.Index(0, size_row) = 2; A.Value(0, size_row) = 2.4;
    DISP(A);

    cout << "Number of entries in the matrix : " << A.GetNonZeros() << endl;

    // if you resize, previous elements are conserved
    A.Resize(m+1,n+1); DISP(A);

    // you can read/write the matrix
    A.Write("mat_binary.dat");
    A.WriteText("mat_ascii.dat");
    A.Clear();
    A.Read("mat_binary.dat"); DISP(A);
    A.Clear();
    A.ReadText("mat_ascii.dat"); DISP(A);

    // you can retrieve pointer to row numbers with GetIndex
    // and values with GetData
    int* row_ptr = A.GetIndex(2);
    cout << "First index " << row_ptr[0] << endl;
    double* value_ptr = A.GetData(2);
    cout << "First value " << value_ptr[0] << endl;
    // you can fill a sparse vectors with pointer, then you call Nullify
    Vector<double, VectSparse> row_vec;
    row_vec.SetData(A.GetRowSize(2), value_ptr, row_ptr);
    // Nullify is necessary to avoid duplicate release of memory
    A.Nullify(2);

    // you can initialize A with null values
    A.Zero(); DISP(A);
    // or with 0, 1, ...
    A.Fill(); DISP(A);
    // or randomly
    A.FillRand(); DISP(A);
    // fill with a given value
    A.Fill(0.01);
    A.Get(2, 3) = 0.5;
    A.Get(0, 1) = -0.6;
    // remove small values (below 0.1 here)
    A.RemoveSmallEntry(0.1); DISP(A);
    // you can set A = I
    A.SetIdentity(); A.Print();

    // with clear, A is an empty matrix with 0 row, 0 column
    A.Clear(); DISP(A);
  }

  {
    // We are testing symmetric matrix.
    int m = 5, n = 5;
    Matrix<double, Symmetric, ArrayRowSymSparse> A(m,n);

    cout << "Number of rows in the matrix " << A.GetM() << endl;
    cout << "Number of columns in the matrix " << A.GetN() << endl;

    // We can add elements by using the method index and value:
    int irow = 3, icol, nb_values_row = 2;
    A.ReallocateRow(irow, nb_values_row);
    // Be careful, only superior part of the matrix is stored (i <= j).
    A.Index(irow, 0) = 4; A.Value(irow, 0) = -2.0;
    A.Index(irow, 1) = 3; A.Value(irow, 1) = 3.0;
    // Values not sorted you have to call.
    A.AssembleRow(irow); DISP(A);

    // with reallocate, A is empty with m columns, n rows
    A.Reallocate(m,n); DISP(A);

    // Be careful, only upper part of the matrix is stored (i <= j).
    irow = 2;
    A.ReallocateRow(irow, nb_values_row);
    A.Index(irow, 0) = 2; A.Value(irow, 0) = 1.0;
    A.Index(irow, 1) = 4; A.Value(irow, 1) = 1.0;
    // Values already sorted, no need to call AssembleRow.
    DISP(A);


    // We can add elements by using the method AddInteraction.
    double value = 2.5; irow = 1; icol = 3;
    A.AddInteraction(irow, icol, value);
    DISP(A);

    // AddInteractionsRow (no need to sort numbers)
    // Elements placed in lower part of the matrix are not added.
    IVect columns(nb_values_row); Vector<double> values(nb_values_row);
    columns(0) = 0; columns(1) = 2;
    values(0) = -3.0; values(1) = 4.6;
    A.AddInteractionRow(irow, nb_values_row, columns, values);
    DISP(A);

    // AddInteractionsColumn (no need to sort numbers)
    // Elements placed in lower part of the matrix are not added.
    int nb_values_col = 3; icol = 0;
    IVect rows(nb_values_col); values.Reallocate(nb_values_col);
    rows(0) = 3; rows(1) = 0; rows(2) = 1;
    values(0) = -1.0; values(1) = 1.5; values(2) = 0.5;
    A.AddInteractionColumn(icol, nb_values_col, rows, values);
    DISP(A);

    // You can directly access A(i,j) and modify it.
    DISP(A(4,2));
    A.Get(4,2) = 2.5; DISP(A(4,2));

    // You can clear a row for symmetric matrix, only upper part of matrix is
    // modified,
    A.ClearRow(3);
    // and resize a row and add new coefficients at the end.
    int size_row = A.GetRowSize(1);
    A.ResizeRow(1, size_row+1);
    A.Index(1, size_row) = 4; A.Value(1, size_row) = 2.4;
    DISP(A);

    // For symmetric matrix, we count non-zero entries only in upper part.
    cout << "Number of entries in the matrix : " << A.GetNonZeros() << endl;

    // if you resize, previous elements are conserved
    A.Resize(m+1,n+1); DISP(A);

    // you can read/write the matrix
    A.Write("mat_binary.dat");
    A.WriteText("mat_ascii.dat");
    A.Clear();
    A.Read("mat_binary.dat"); DISP(A);
    A.Clear();
    A.ReadText("mat_ascii.dat"); DISP(A);

    // you can retrieve pointer to row numbers with GetIndex
    // and values with GetData
    int* row_ptr = A.GetIndex(2);
    cout << "First index " << row_ptr[0] << endl;
    double* value_ptr = A.GetData(2);
    cout << "First value " << value_ptr[0] << endl;
    // you can fill a sparse vectors with pointer, then you call Nullify
    Vector<double, VectSparse> row_vec;
    row_vec.SetData(A.GetRowSize(2), value_ptr, row_ptr);
    // Nullify is necessary to avoid duplicate release of memory
    A.Nullify(2);

    // you can initialize A with null values
    A.Zero(); DISP(A);
    // or with 0, 1, ...
    A.Fill(); DISP(A);
    // or randomly
    A.FillRand(); DISP(A);
    // fill with a given value
    A.Fill(0.01);
    A.Get(2, 3) = 0.5;
    A.Get(0, 1) = -0.6;
    // remove small values (below 0.1 here)
    A.RemoveSmallEntry(0.1); DISP(A);
    // you can set A = I
    A.SetIdentity(); A.Print();

    // with clear, A is an empty matrix with 0 row, 0 column
    A.Clear(); DISP(A);
  }

  return 0;
}
