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


#ifndef SELDON_FILE_ERRORS_CXX

#include "Errors.hxx"

namespace Seldon
{


  ///////////
  // ERROR //
  ///////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  Error::Error(string function, string comment):
    description_("ERROR!\nAn undefined error occurred"),
    function_(function), comment_(comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }


  //! Alternative constructor.
  /*! Error associated with a description, a function and a comment.
    \param description short description of the error.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  Error::Error(string description, string function, string comment):
    description_("ERROR!\n" + description),
    function_(function), comment_(comment)
  {
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  /*!
    \note Empty.
  */
  Error::~Error()
  {
  }


  /***********
   * METHODS *
   ***********/


  //! Delivers information about the error.
  /*! Displays available information, i.e. the error description, the function
    and/or a comment.
  */
  string Error::What()
  {
    string message(description_);
    if (!function_.empty())
      message += " in " + function_;
    message += ".\n";
    if (!comment_.empty())
      message += "   " + comment_;
    return message;
  }


  //! Delivers information about the error.
  /*! Displays available information, i.e. the error description, the function
    and/or a comment.
  */
  void Error::CoutWhat()
  {
    cout << this->What() << endl;
  }



  ///////////////
  // UNDEFINED //
  ///////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param[in] function function with which the error is associated.
    \param[in] comment comment associated with the error.
  */
  Undefined::Undefined(string function, string comment):
  Error("", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }


  //! Delivers information about the error.
  /*! Displays available information, i.e.
    the error description, the function and/or the comment.
  */
  string Undefined::What()
  {
    string message;
    if (!this->function_.empty())
      message = this->function_;
    else
      message = "A function or a method";
    message += " is undefined.\nEither its implementation is missing,"
      + string(" or it does not make sense or it is impossible ")
      + "to implement it.\n";
    if (!this->comment_.empty())
      message += "   " + this->comment_;
    return message;
  }


  ///////////////////
  // WRONGARGUMENT //
  ///////////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  WrongArgument::WrongArgument(string function, string comment):
  Error("Wrong argument given to ", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }


  //! Delivers information about the error.
  /*! Displays available information, i.e.
    the error description, the function and/or the comment.
  */
  string WrongArgument::What()
  {
    string message(this->description_);
    if (!this->function_.empty())
      message += this->function_;
    message += ".\n";
    if (!this->comment_.empty())
      message += "   " + this->comment_;
    return message;
  }


  //////////////
  // NOMEMORY //
  //////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  NoMemory::NoMemory(string function, string comment):
    Error("Out of memory", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }



  //////////////
  // WRONGDIM //
  //////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  WrongDim::WrongDim(string function, string comment):
    Error("Wrong dimensions involved", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }



  ////////////////
  // WRONGINDEX //
  ////////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  WrongIndex::WrongIndex(string function, string comment):
    Error("Index out of range", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }



  //////////////
  // WRONGROW //
  //////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  WrongRow::WrongRow(string function, string comment):
    Error("Row index out of range", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }



  //////////////
  // WRONGCOL //
  //////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  WrongCol::WrongCol(string function, string comment):
    Error("Column index out of range", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }



  /////////////
  // IOERROR //
  /////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  IOError::IOError(string function, string comment):
    Error("Error while performing an I/O operation", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }



  /////////////////
  // LAPACKERROR //
  /////////////////


  //! Main constructor.
  /*! Error associated with a diagnostic integer, a function and a comment.
    \param info Lapack diagnostic integer.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  LapackError::LapackError(int info, string function, string comment):
  Error("Error returned by Lapack", function, comment), info_(info)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }


  //! Delivers information about the error.
  /*! Displays available information, i.e.
    the error description, the function and/or the comment.
  */
  string LapackError::What()
  {
    string message(this->description_);
    if (!this->function_.empty())
      message += " in " + this->function_;
    message += ".\n";
    if (!this->comment_.empty())
      message += "   " + this->comment_;
    message += "   Diagnostic integer (\"info\"): " + to_str(info_)
      + ".\n";
    return message;
  }


  ////////////////
  // LAPACKINFO //
  ////////////////


  LapackInfo::LapackInfo(int info): info_(info)
  {
  }


  LapackInfo::operator int ()
  {
    return info_;
  }


  int LapackInfo::GetInfo()
  {
    return info_;
  }


  int& LapackInfo::GetInfoRef()
  {
    return info_;
  }


#ifndef SELDON_WITH_COMPILED_LIBRARY
  LapackInfo lapack_info(0);
#endif


} // namespace Seldon.

#define SELDON_FILE_ERRORS_CXX
#endif
