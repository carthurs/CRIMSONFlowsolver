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


#ifndef VERDANDI_FILE_SHARE_ERROR_HXX


namespace Verdandi
{


    /*! \brief This class serves for exceptions raised when an error is
      found. */
    /*! All Verdandi exceptions are supposed to inherit from this class. */
    class Error
    {
    protected:
        //! Message describing the exception type.
        string description_;
        //! Name of the function in which the error occurred.
        string function_;
        //! A comment about the error.
        string comment_;

    public:
        // Constructors.
        Error(string function, string comment) throw();
        Error(string description, string function, string comment) throw();

        // Destructor.
        virtual ~Error() throw();

        virtual  string GetName() const;
        virtual string What();
        void CoutWhat();
    };


    /*! \brief This class serves for exceptions raised when an error is found
      in a configuration. */
    class ErrorConfiguration: public Error
    {
    public:
        // Constructor.
        ErrorConfiguration(string function, string comment) throw();

        string GetName() const;
    };


    /*! \brief This class serves for exceptions raised when an input/output
      error occurs. */
    class ErrorIO: public Error
    {
    public:
        // Constructor.
        ErrorIO(string function, string comment) throw();

        string GetName() const;
    };


    /*! \brief This class serves for exceptions raised when an error occurs
      during some data processing. */
    class ErrorProcessing: public Error
    {
    public:
        // Constructor.
        ErrorProcessing(string function, string comment) throw();

        string GetName() const;
    };


    /*! \brief This class serves for exceptions raised when an undefined
      function or method is called. */
    class ErrorUndefined: public Error
    {
    public:
        // Constructor.
        ErrorUndefined(string function, string comment) throw();

        virtual string What();
        string GetName() const;
    };


    /*! \brief This class serves for exceptions raised when a function or
      a method is called with an erroneous argument. */
    class ErrorArgument: public Error
    {
    public:
        // Constructor.
        ErrorArgument(string function, string comment) throw();

        string GetName() const;
    };


    /*! \brief This class serves for exceptions raised when an undefined
      function or method in a Python module is called. */
    class ErrorPythonUndefined: public Error
    {
    protected:
        //! Name of the Python function.
        string function_name_;
        //! List of arguments of the function.
        string arguments_;
        //! Name of the Python module.
        string module_;
    public:
        // Constructor.
        ErrorPythonUndefined(string function, string function_name,
                             string arguments,
                             string module,
                             string comment) throw();

        // Destructor.
        ~ErrorPythonUndefined() throw ();

        virtual string What();
        string GetName() const;
    };


} // namespace Verdandi.


#define VERDANDI_FILE_SHARE_ERROR_HXX
#endif
