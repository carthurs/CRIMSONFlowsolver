#define SELDON_DEBUG_LEVEL_3

#include "Seldon.hxx"
using namespace Seldon;

int main()
{

  TRY;

  Matrix<double> A(3, 3);
  Vector<double> U, V;

  A.SetIdentity();
  A(0, 1) = -1.0;

  U.Reallocate(A.GetN());
  U.Fill();

  U.Print();

  V.Reallocate(A.GetM());
  Mlt(2.0, A, U, V);

  cout << V << endl;

  END;

  return 0;

}
