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


#define SELDON_WITH_BLAS
#include "Seldon.hxx"
using namespace Seldon;

int main()
{

  cout << "Seldon: compilation test with Blas" << endl;

  Vector<double> U(3), V(3);
  U.Fill(1.3);
  V.Fill();

  cout << "U = " << U << endl;
  cout << "V = " << V << endl;

  // 2.0 * U + V -> V
  Add(2.0, U, V);

  // Note: 'Add' calls the Blas function 'daxpy'.

  cout << "2. * U + V = " << V << endl;

  return 0;

}
