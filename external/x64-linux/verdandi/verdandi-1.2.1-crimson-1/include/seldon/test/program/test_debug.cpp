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


#define SELDON_DEBUG_LEVEL_4
#define SELDON_DEFAULT_ALLOCATOR NewAlloc
#define SELDON_WITH_ABORT

#include "Seldon.hxx"

using namespace Seldon;

// testing class Vector without Blas

int main()
{
  Vector<double> V(0);
  // Vector<double> U(-2);
  V.Reallocate(0);
  V.Reallocate(-2);
  V(0) = 1.0;
  V.Clear();
}
