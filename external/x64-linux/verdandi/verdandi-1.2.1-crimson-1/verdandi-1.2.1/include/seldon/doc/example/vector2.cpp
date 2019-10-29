#define SELDON_DEBUG_LEVEL_4
#include "Seldon.hxx"
using namespace Seldon;

#include "vector/Vector2.cxx"

int main()
{

  TRY;

  Vector<int> length(3);
  length(0) = 2;
  length(1) = 3;
  length(2) = 7;
  Vector2<double> V(length);

  // Fills all inner vectors with 2.
  V.Fill(2.);
  // Access to the second inner vector, to fill it with 5.
  V(1).Fill(5.);

  V.Print();
  cout << "First element of the second inner vector: " << V(1, 0) << endl;
  // Note that V(1)(0) would have returned the same element.

  Vector<double> inner_vector(4);
  inner_vector.Fill();
  // Appends a new inner vector.
  V.PushBack(inner_vector);
  V.Print();

  cout << "After setting to -10 the second element of the last inner vector:"
       << endl;
  V(3, 1) = -10.;
  V.Print();

  END;

  return 0;

}
