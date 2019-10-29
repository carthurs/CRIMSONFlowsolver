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


#ifndef OPS_FILE_ERROR_CXX


#include "Error.hxx"


namespace Ops
{


  //! Constructor.
  /*! Error associated with both a function and a comment.
    \param[in] function function in which the error occurred.
    \param[in] comment comment associated with the error.
  */
  Error::Error(string function = "", string comment = ""):
    function_(function), comment_(comment)
  {
#ifdef OPS_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }


  //! Destructor.
  Error::~Error()
  {
  }


  //! Delivers information about the error.
  /*! \return The available information: the function and/or the comment.
   */
  string Error::What()
  {
    string message;
    if (!function_.empty())
      message += "Error in Ops::" + function_ + ".\n";
    if (!comment_.empty())
      message += "   " + comment_;
    return message;
  }


  //! Delivers information about the error.
  /*! Displays available information: the function and/or the comment.
   */
  void Error::CoutWhat()
  {
    std::cout << this->What() << std::endl;
  }


} // namespace Ops.


#define OPS_FILE_ERROR_CXX
#endif
