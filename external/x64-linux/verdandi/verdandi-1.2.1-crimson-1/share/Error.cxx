// Copyright (C) 2008-2009, INRIA
// Author(s): Vivien Mallet
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


#ifndef VERDANDI_FILE_SHARE_ERROR_CXX


#include "Error.hxx"


namespace Verdandi
{


    ///////////
    // ERROR //
    ///////////


    /****************
     * CONSTRUCTORS *
     ****************/


    //! Main constructor.
    /*! Error associated with both a function and a comment.
      \param[in] function function in which the error occurred.
      \param[in] comment comment associated with the error.
    */
    Error::Error(string function = "", string comment = "") throw():
        description_("An undefined error occurred"),
        function_(function), comment_(comment)
    {
#ifdef VERDANDI_WITH_ABORT
        this->CoutWhat();
        Logger::SetStdout(false);
        Logger::Log<VERDANDI_EXCEPTION_LOGGING_LEVEL>(*this, this->What());
        abort();
#else
        Logger::Log<VERDANDI_EXCEPTION_LOGGING_LEVEL>(*this, this->What());
#endif
    }


    //! Alternative constructor.
    /*! Error associated with a description, a function and a comment.
      \param[in] description short description of the error.
      \param[in] function function in which the error occurred.
      \param[in] comment comment associated with the error.
    */
    Error::Error(string description, string function, string comment) throw():
        description_(description), function_(function), comment_(comment)
    {
    }


    /**************
     * DESTRUCTOR *
     **************/


    //! Destructor.
    /*!
      \note Empty.
    */
    Error::~Error() throw()
    {
    }


    /***********
     * METHODS *
     ***********/


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    string Error::GetName() const
    {
        return "Error";
    }


    //! Delivers information about the error.
    /*! Displays available information, i.e. the error description, the
      function and/or the comment.
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
    /*! Displays available information, i.e. the error description, the
      function and/or the comment.
    */
    void Error::CoutWhat()
    {
        cout << this->What() << endl;
    }


    ////////////////////////
    // ERRORCONFIGURATION //
    ////////////////////////


    //! Main constructor.
    /*! Error associated with both a function and a comment.
      \param[in] function function in which the error occurred.
      \param[in] comment comment associated with the error.
    */
    ErrorConfiguration::ErrorConfiguration(string function = "",
                                           string comment = "")
        throw():
        Error("Error while reading a configuration file", function, comment)
    {
#ifdef VERDANDI_WITH_ABORT
        this->CoutWhat();
        Logger::SetStdout(false);
        Logger::Log<VERDANDI_EXCEPTION_LOGGING_LEVEL>(*this, this->What());
        abort();
#else
        Logger::Log<VERDANDI_EXCEPTION_LOGGING_LEVEL>(*this, this->What());
#endif
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    string ErrorConfiguration::GetName() const
    {
        return "ErrorConfiguration";
    }


    /////////////
    // ERRORIO //
    /////////////


    //! Main constructor.
    /*! Error associated with both a function and a comment.
      \param[in] function function in which the error occurred.
      \param[in] comment comment associated with the error.
    */
    ErrorIO::ErrorIO(string function = "", string comment = "")
        throw():
        Error("Error while performing a I/O operation", function, comment)
    {
#ifdef VERDANDI_WITH_ABORT
        this->CoutWhat();
        Logger::SetStdout(false);
        Logger::Log<VERDANDI_EXCEPTION_LOGGING_LEVEL>(*this, this->What());
        abort();
#else
        Logger::Log<VERDANDI_EXCEPTION_LOGGING_LEVEL>(*this, this->What());
#endif
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    string ErrorIO::GetName() const
    {
        return "ErrorIO";
    }


    /////////////////////
    // ERRORPROCESSING //
    /////////////////////


    //! Main constructor.
    /*! Error associated with both a function and a comment.
      \param[in] function function in which the error occurred.
      \param[in] comment comment associated with the error.
    */
    ErrorProcessing
    ::ErrorProcessing(string function = "", string comment = "") throw():
        Error("Error while processing data", function, comment)
    {
#ifdef VERDANDI_WITH_ABORT
        this->CoutWhat();
        Logger::SetStdout(false);
        Logger::Log<VERDANDI_EXCEPTION_LOGGING_LEVEL>(*this, this->What());
        abort();
#else
        Logger::Log<VERDANDI_EXCEPTION_LOGGING_LEVEL>(*this, this->What());
#endif
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    string ErrorProcessing::GetName() const
    {
        return "ErrorProcessing";
    }


    ////////////////////
    // ERRORUNDEFINED //
    ////////////////////


    //! Main constructor.
    /*! Error associated with both a function and a comment.
      \param[in] function function in which the error occurred.
      \param[in] comment comment associated with the error.
    */
    ErrorUndefined::ErrorUndefined(string function = "", string comment = "")
        throw():
        Error("Undefined function", function, comment)
    {
#ifdef VERDANDI_WITH_ABORT
        this->CoutWhat();
        Logger::SetStdout(false);
        Logger::Log<VERDANDI_EXCEPTION_LOGGING_LEVEL>(*this, this->What());
        abort();
#else
        Logger::Log<VERDANDI_EXCEPTION_LOGGING_LEVEL>(*this, this->What());
#endif
    }


    //! Delivers information about the error.
    /*! Displays available information, i.e. the error description, the
      function and/or the comment.
    */
    string ErrorUndefined::What()
    {
        string message(description_);
        if (!function_.empty())
            message += ": " + function_;
        message += ".\n";
        if (!comment_.empty())
            message += "   " + comment_;
        return message;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    string ErrorUndefined::GetName() const
    {
        return "ErrorUndefined";
    }


    ///////////////////
    // ERRORARGUMENT //
    ///////////////////


    //! Main constructor.
    /*! Error associated with both a function and a comment.
      \param[in] function function in which the error occurred.
      \param[in] comment comment associated with the error.
    */
    ErrorArgument::ErrorArgument(string function = "", string comment = "")
        throw():
        Error("Wrong arguments", function, comment)
    {
#ifdef VERDANDI_WITH_ABORT
        this->CoutWhat();
        Logger::SetStdout(false);
        Logger::Log<VERDANDI_EXCEPTION_LOGGING_LEVEL>(*this, this->What());
        abort();
#else
        Logger::Log<VERDANDI_EXCEPTION_LOGGING_LEVEL>(*this, this->What());
#endif
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    string ErrorArgument::GetName() const
    {
        return "ErrorArgument";
    }


    //////////////////////////
    // ERRORPYTHONUNDEFINED //
    //////////////////////////


    //! Main constructor.
    /*! Error associated with both a function and a comment.
      \param[in] function function in which the error occurred.
      \param[in] function_name name of the python function called.
      \param[in] arguments arguments that must be passed to the function.
      \param[in] module Python module where the function is searched.
      \param[in] comment comment associated with the error.
    */
    ErrorPythonUndefined::ErrorPythonUndefined(string function = "",
                                               string function_name = "",
                                               string arguments = "",
                                               string module = "",
                                               string comment = "")
        throw():
        Error("Call to an undefined Python function", function, comment)
    {
        function_name_ = function_name;
        arguments_ = arguments;
        module_ = module;

#ifdef VERDANDI_WITH_ABORT
        this->CoutWhat();
        Logger::SetStdout(false);
        Logger::Log<VERDANDI_EXCEPTION_LOGGING_LEVEL>(*this, this->What());
        abort();
#else
        Logger::Log<VERDANDI_EXCEPTION_LOGGING_LEVEL>(*this, this->What());
#endif
    }


    //! Destructor.
    /*!
      \note Empty.
    */
    ErrorPythonUndefined::~ErrorPythonUndefined() throw()
    {
    }

    //! Delivers information about the error.
    /*! Displays available information, i.e. the error description, the
      function and/or the comment.
    */
    string ErrorPythonUndefined::What()
    {
        string message(description_);
        if (!function_.empty())
            message += " in " + function_;
        message += ".\n";
        message += "   The Python function '" + function_name_;
        message += arguments_;
        message += "' is either not defined in module \"" + module_;
        message += "\" or has wrong arguments.";
        if (!comment_.empty())
            message += "\n    " + comment_;
        return message;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    string ErrorPythonUndefined::GetName() const
    {
        return "ErrorPythonUndefined";
    }


} // namespace Verdandi.


#define VERDANDI_FILE_SHARE_ERROR_CXX
#endif

