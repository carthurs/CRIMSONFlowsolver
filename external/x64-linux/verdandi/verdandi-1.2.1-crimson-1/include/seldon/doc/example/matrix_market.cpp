#define SELDON_DEBUG_LEVEL_2

#include "Seldon.hxx"
using namespace Seldon;

#include "matrix_sparse/IOMatrixMarket.cxx"

int main(int argc, char** argv)
{

  TRY;

  Matrix<double, General, ColSparse> A;

  ReadHarwellBoeing(argv[1], A);

  A.Print();

  WriteHarwellBoeing(A, "result.rua");

  END;

  return 0;

}
