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


#ifndef VERDANDI_FILE_OUTPUTSAVER_VARIABLE_HXX


namespace Verdandi
{


    //////////////
    // VARIABLE //
    /////////////


    //! This class enables to define a variable for the output saver.
    class Variable
    {

    private:

        /*** Main components ***/

        //! Saving format (e.g., "binary", "text").
        string mode_;
        //! Path to the output file.
        string file_;
        //! Boolean to indicate if the file has to be emptied.
        bool has_to_empty_file_;

    public:

        /*** Constructors and destructor ***/

        Variable();
        Variable(string mode, string file, bool has_to_empty_file = false);
        Variable(const Variable& variable);
        ~Variable();

        /*** Access methods ***/

        void SetMode(string mode);
        void SetFile(string file);
        void HasToEmptyFile(bool has_to_empty_file);

        string GetMode() const;
        string GetFile() const;
        bool  HasToEmptyFile() const;

        void Display() const;

    };


} // namespace Verdandi.


#define VERDANDI_FILE_OUTPUTSAVER_VARIABLE_HXX
#endif
