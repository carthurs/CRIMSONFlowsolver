#define SELDON_DEBUG_LEVEL_4

#include "Seldon.hxx"
using namespace Seldon;

int main()
{

  TRY;

  Matrix<double> A(3, 3);

  A.Zero();
  A(0, 3) = 2.0;

  END;

  cout << "The program should not reach this point..." << endl;

  return 0;

}
