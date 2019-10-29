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


#ifndef SELDON_FILE_STORAGE_CXX

#include "Storage.hxx"

namespace Seldon
{


  //////////////////////
  // GENERAL MATRICES //
  //////////////////////


  inline int ColMajor::GetFirst(int i, int j)
  {
    return j;
  }
  inline int ColMajor::GetSecond(int i, int j)
  {
    return i;
  }


  inline int RowMajor::GetFirst(int i, int j)
  {
    return i;
  }
  inline int RowMajor::GetSecond(int i, int j)
  {
    return j;
  }



  /////////////
  // VECTORS //
  /////////////


  class VectFull
  {
  };


  class VectSparse
  {
  };


  class Collection
  {
  };


  class DenseSparseCollection
  {
  };


  class PETScSeq
  {
  };


  class PETScPar
  {
  };


  class PETScSeqDense
  {
  };


  class PETScMPIDense
  {
  };


  class PETScMPIAIJ
  {
  };


  ////////////
  // SPARSE //
  ////////////


  inline int ColSparse::GetFirst(int i, int j)
  {
    return j;
  }
  inline int ColSparse::GetSecond(int i, int j)
  {
    return i;
  }


  inline int RowSparse::GetFirst(int i, int j)
  {
    return i;
  }
  inline int RowSparse::GetSecond(int i, int j)
  {
    return j;
  }


  inline int ColComplexSparse::GetFirst(int i, int j)
  {
    return j;
  }
  inline int ColComplexSparse::GetSecond(int i, int j)
  {
    return i;
  }


  inline int RowComplexSparse::GetFirst(int i, int j)
  {
    return i;
  }
  inline int RowComplexSparse::GetSecond(int i, int j)
  {
    return j;
  }


  inline int ColSymSparse::GetFirst(int i, int j)
  {
    return j;
  }
  inline int ColSymSparse::GetSecond(int i, int j)
  {
    return i;
  }


  inline int RowSymSparse::GetFirst(int i, int j)
  {
    return i;
  }
  inline int RowSymSparse::GetSecond(int i, int j)
  {
    return j;
  }


  inline int ColSymComplexSparse::GetFirst(int i, int j)
  {
    return j;
  }
  inline int ColSymComplexSparse::GetSecond(int i, int j)
  {
    return i;
  }


  inline int RowSymComplexSparse::GetFirst(int i, int j)
  {
    return i;
  }
  inline int RowSymComplexSparse::GetSecond(int i, int j)
  {
    return j;
  }



  ///////////////
  // SYMMETRIC //
  ///////////////


  inline int ColSymPacked::GetFirst(int i, int j)
  {
    return j;
  }
  inline int ColSymPacked::GetSecond(int i, int j)
  {
    return i;
  }


  inline int RowSymPacked::GetFirst(int i, int j)
  {
    return i;
  }
  inline int RowSymPacked::GetSecond(int i, int j)
  {
    return j;
  }


  inline int ColSym::GetFirst(int i, int j)
  {
    return j;
  }
  inline int ColSym::GetSecond(int i, int j)
  {
    return i;
  }


  inline int RowSym::GetFirst(int i, int j)
  {
    return i;
  }
  inline int RowSym::GetSecond(int i, int j)
  {
    return j;
  }



  ///////////////
  // HERMITIAN //
  ///////////////


  inline int ColHerm::GetFirst(int i, int j)
  {
    return j;
  }
  inline int ColHerm::GetSecond(int i, int j)
  {
    return i;
  }


  inline int RowHerm::GetFirst(int i, int j)
  {
    return i;
  }
  inline int RowHerm::GetSecond(int i, int j)
  {
    return j;
  }


  inline int ColHermPacked::GetFirst(int i, int j)
  {
    return j;
  }
  inline int ColHermPacked::GetSecond(int i, int j)
  {
    return i;
  }


  inline int RowHermPacked::GetFirst(int i, int j)
  {
    return i;
  }
  inline int RowHermPacked::GetSecond(int i, int j)
  {
    return j;
  }



  ////////////////
  // TRIANGULAR //
  ////////////////


  inline int ColUpTriang::GetFirst(int i, int j)
  {
    return j;
  }
  inline int ColUpTriang::GetSecond(int i, int j)
  {
    return i;
  }
  inline bool ColUpTriang::UpLo()
  {
    return true;
  }


  inline int ColLoTriang::GetFirst(int i, int j)
  {
    return j;
  }
  inline int ColLoTriang::GetSecond(int i, int j)
  {
    return i;
  }
  inline bool ColLoTriang::UpLo()
  {
    return false;
  }


  inline int RowUpTriang::GetFirst(int i, int j)
  {
    return i;
  }
  inline int RowUpTriang::GetSecond(int i, int j)
  {
    return j;
  }
  inline bool RowUpTriang::UpLo()
  {
    return true;
  }


  inline int RowLoTriang::GetFirst(int i, int j)
  {
    return i;
  }
  inline int RowLoTriang::GetSecond(int i, int j)
  {
    return j;
  }
  inline bool RowLoTriang::UpLo()
  {
    return false;
  }


  inline int ColUpTriangPacked::GetFirst(int i, int j)
  {
    return j;
  }
  inline int ColUpTriangPacked::GetSecond(int i, int j)
  {
    return i;
  }
  inline bool ColUpTriangPacked::UpLo()
  {
    return true;
  }


  inline int ColLoTriangPacked::GetFirst(int i, int j)
  {
    return j;
  }
  inline int ColLoTriangPacked::GetSecond(int i, int j)
  {
    return i;
  }
  inline bool ColLoTriangPacked::UpLo()
  {
    return false;
  }


  inline int RowUpTriangPacked::GetFirst(int i, int j)
  {
    return i;
  }
  inline int RowUpTriangPacked::GetSecond(int i, int j)
  {
    return j;
  }
  inline bool RowUpTriangPacked::UpLo()
  {
    return true;
  }


  inline int RowLoTriangPacked::GetFirst(int i, int j)
  {
    return i;
  }
  inline int RowLoTriangPacked::GetSecond(int i, int j)
  {
    return j;
  }
  inline bool RowLoTriangPacked::UpLo()
  {
    return false;
  }


} // namespace Seldon.

#define SELDON_FILE_STORAGE_CXX
#endif
