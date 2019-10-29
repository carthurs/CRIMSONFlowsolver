// Copyright (C) 2011 Marc Duruflé
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


#ifndef SELDON_FILE_SOLVER_ORDERING_HXX
#define SELDON_FILE_SOLVER_ORDERING_HXX


namespace Seldon
{

  //! Basic class grouping different ordering strategies.
  class SparseMatrixOrdering
  {
  public :
    // Supported orderings.
    enum {IDENTITY, REVERSE_CUTHILL_MCKEE, PORD,
	  SCOTCH, METIS, AMD, COLAMD, QAMD, USER, AUTO};
  };


  template<class T, class Prop, class Storage, class Allocator,
	   class Tint, class Alloc>
  void FindSparseOrdering(Matrix<T, Prop, Storage, Allocator>& A,
			  Vector<Tint, VectFull, Alloc>& num, int type);

} // namespace Seldon.


#endif
