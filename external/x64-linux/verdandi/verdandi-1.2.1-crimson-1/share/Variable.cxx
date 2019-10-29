// Copyright (C) 2008-2009 INRIA
// Author(s): Marc Fragu
//
// This file is part of the data assimilation library Verdandi.
//
// Verdandi is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Verdandi. If not, see http://www.gnu.org/licenses/.
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_FILE_OUTPUTSAVER_VARIABLE_CXX


#include "Variable.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    Variable::Variable()
    {
    }


    //! Main constructor.
    /*!
      \param[in] mode saving format (e.g., "binary", "text").
      \param[in] file path to the output file.
      \param[in] has_to_empty_file Boolean that indicates if the file has to
      be emptied.
    */
    Variable::Variable(string mode, string file, bool has_to_empty_file):
        mode_(mode), file_(file), has_to_empty_file_(has_to_empty_file)
    {
    }


    //! Copy constructor.
    /*!
      \param[in] variable instance to be copied.
    */
    Variable::Variable(const Variable& variable)
    {
        mode_ = variable.GetMode();
        file_ = variable.GetFile();
    }


    //! Destructor.
    Variable::~Variable()
    {
    }


    ///////////////////
    // ACCESS METHOD //
    ///////////////////


    //! Mode accessor.
    /*! Sets the saving mode.
      \param[in] mode saving format (e.g., "binary", "text").
    */
    void Variable::SetMode (string mode)
    {
        mode_ = mode;
    }


    //! Filename accessor.
    /*! Sets the path to the output file.
      \param[in] file path to the output file.
    */
    void Variable::SetFile(string file)
    {
        file_ = file;
    }


    //! Empty accessor.
    /*! Sets the value of the boolean empty.
      \param[in] has_to_empty_file Boolean that indicates if the file has to
      be emptied.
    */
    void Variable::HasToEmptyFile(bool has_to_empty_file)
    {
        has_to_empty_file_ = has_to_empty_file;
    }


    //! Mode accessor.
    /*! Returns the saving mode.
      \return The saving format.
    */
    string Variable::GetMode() const
    {
        return mode_;
    }


    //! File accessor.
    /*! Returns the path to the output file.
      \return The path to the output file.
    */
    string Variable::GetFile() const
    {
        return file_;
    }


    //! Empty accessor.
    /*! Indicates if the file has to be emptied.
      \return Boolean that indicates if the file has to be emptied.
    */
    bool Variable::HasToEmptyFile() const
    {
        return has_to_empty_file_;
    }


    //! Displays on screen the saving mode and the output filename.
    void Variable::Display() const
    {
        cout << "Mode: " << mode_ << "\n";
        cout << "File: " << file_ << endl;
    }


}


#define VERDANDI_FILE_OUTPUTSAVER_VARIABLE_CXX
#endif

