// Copyright (C) 2008-2009 INRIA
// Author(s): Marc Fragu, Vivien Mallet
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


#ifndef VERDANDI_FILE_OUTPUTSAVER_OUTPUTSAVER_HXX


#include "Variable.hxx"

#include <iostream>
#include <fstream>
#include <string>


namespace Verdandi
{


    /////////////////
    // OUTPUTSAVER //
    /////////////////



    //! This class provides convenient methods to save variables on disk.
    class OutputSaver
    {

    private:

        //! Time interval between two saves.
        double save_period_;
        //! Tolerance on the time interval.
        double time_tolerance_;

        //! Default mode.
        string mode_;

        //! Default mode for scalar variables.
        string mode_scalar_;

        //! Stores the variables properties.
        map<string, Variable> variable_list_;

        //! Boolean to indicate if the output saver is active or not.
        bool is_active_;

    public:

        /*** Constructors and destructor ***/

        OutputSaver();
        OutputSaver(string configuration_file, string method_name);
        void Initialize(string configuration_file, string method_name);
        void Initialize(VerdandiOps& configuration);
        void Activate();
        void Deactivate();

        ~OutputSaver();

        string GetName() const;

        /*** Methods ***/

        template <class S>
        void Save(const S& x, double time, string variable_name);

        template <class S>
        void Save(const S& x, string variable_name);
#ifdef VERDANDI_WITH_PETSC
        template <class T, class Allocator>
        void Save(const Vector<T, PETScPar, Allocator>& x,
                  string variable_name);
#endif

        template <class S>
        void WriteText(const S& x, string file_name) const;
        template <class S>
        void WriteBinary(const S& x, string file_name) const;
        template <class T, class Prop, class Allocator>
        void WriteBinary(const Matrix<T, Prop, RowSparse, Allocator>& x,
                         string file_name) const;

        template <class T, class Prop, class Allocator>
        void WriteBinary(const Matrix<T, Prop, ColSparse, Allocator>& x,
                         string file_name) const;


        template <class S>
        void Empty(string variable_name);
        void Empty(string variable_name);
        void Empty();

        bool IsVariable(string variable_name) const;
        void DisplayVariableList() const;

    private:

        void SetVariable(VerdandiOps& configuration,
                         string generic_path,
                         string default_mode,
                         string variable_name);
        template <class S>
        void SetVariable(Variable& variable);
        void SetVariableFile(Variable& variable);

    };


} // namespace Verdandi.


#define VERDANDI_FILE_OUTPUTSAVER_OUTPUTSAVER_HXX
#endif

