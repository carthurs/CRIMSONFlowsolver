// Copyright (C) 2010, Vivien Mallet
//
// This file is part of Ops, a library for parsing Lua configuration files.
//
// Ops is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Ops is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Ops. If not, see http://www.gnu.org/licenses/.


#ifndef OPS_FILE_ERROR_HXX


#ifdef OPS_WITH_ABORT
#include <cstdlib>
#endif


namespace Ops
{


  //! An instance of this class is thrown when an error is detected by Ops.
  class Error
  {
  protected:
    //! Name of the function in which the error occurred.
    string function_;
    //! A comment about the error.
    string comment_;

  public:
    // Constructor.
    Error(string function, string comment);

    // Destructor.
    virtual ~Error();

    virtual string What();
    void CoutWhat();
  };


} // namespace Ops.


#define OPS_FILE_ERROR_HXX
#endif
