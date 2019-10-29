#define SELDON_DEBUG_LEVEL_4
#include "Seldon.hxx"
using namespace Seldon;

int main()
{

  TRY;

  /*** Full matrix ***/

  Matrix<double> A(4, 6);
  A.Fill();
  cout << "Complete matrix:" << endl;
  A.Print();

  Vector<int> row_list(2);
  row_list(0) = 1;
  row_list(1) = 2;
  Vector<int> column_list(3);
  column_list(0) = 0;
  column_list(1) = 1;
  column_list(2) = 5;
  SubMatrix<Matrix<double> > SubA(A, row_list, column_list);

  cout << "Sub-matrix:" << endl;
  SubA.Print();

  // Basic operations are supported, but they are slow (Blas/Lapack will not
  // be called).
  Vector<double> X(3), Y(2);
  X.Fill();
  cout << "Multiplied by X = [" << X << "]:" << endl;
  Mlt(SubA, X, Y);
  Y.Print();

  /*** Symmetric matrix ***/

  Matrix<double, General, ColSymPacked> B(4);
  B.Fill();
  cout << "\nComplete matrix:" << endl;
  B.Print();

  row_list(0) = 1;
  row_list(1) = 3;
  column_list(0) = 0;
  column_list(1) = 2;
  column_list(2) = 3;
  SubMatrix<Matrix<double, General, ColSymPacked> >
    SubB(B, row_list, column_list);

  cout << "Sub-matrix (no more symmetric):" << endl;
  SubB.Print();

  // Assignments in the sub-matrix.
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 3; j++)
      SubB(i, j) = -1.;

  // 'B' will remain symmetric.
  cout << "Complete matrix after the sub-matrix is filled with -1:" << endl;
  B.Print();

  END;

  return 0;

}
