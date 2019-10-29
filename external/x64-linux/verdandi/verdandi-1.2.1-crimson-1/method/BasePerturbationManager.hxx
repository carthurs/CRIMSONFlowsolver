// Copyright (C) 2010 INRIA
// Author(s): Vivien Mallet, Anne Tilloy, Kévin Charpentier
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


#ifndef VERDANDI_FILE_METHOD_BASEPERTURBATIONMANAGER_HXX

#include "seldon/vector/VectorCollection.hxx"

namespace Verdandi
{


    /////////////////////////////
    // BASEPERTURBATIONMANAGER //
    /////////////////////////////


    //! This class generates and applies pertubations.
    template<class Derived>
    class BasePerturbationManager: public VerdandiBase
    {
    public:

        /*** Constructors and destructor ***/

        BasePerturbationManager();
        BasePerturbationManager(string configuration_file);
        ~BasePerturbationManager();

        /*** Methods ***/

        void Initialize(string configuration_file);

        template <class T0, class Prop0, class Allocator0,
                  class T1, class Allocator1>
        void Sample(string pdf,
                    Matrix<T0, Prop0, RowSymPacked, Allocator0>& variance,
                    Vector<double, VectFull>& parameter,
                    Vector<double, VectFull>& correlation,
                    Vector<T1, VectFull, Allocator1>& output);

        template <class T0, class Prop0, class Allocator0,
                  class T1, class Allocator1>
        void Sample(string pdf,
                    Matrix<T0, Prop0, RowSymPacked, Allocator0>& variance,
                    Vector<double, VectFull>& parameter,
                    Vector<double, VectFull>& correlation,
                    Vector<T1, Collection, Allocator1>& output);

        template <class T0,
                  class T1, class Allocator1>
        void Sample(string pdf, T0 variance,
                    Vector<double, VectFull>& parameter,
                    Vector<double, VectFull>& correlation,
                    Vector<T1, VectFull, Allocator1>& output);

        template <class T0,
                  class T1, class Allocator1>
        void Sample(string pdf, T0 variance,
                    Vector<double, VectFull>& parameter,
                    Vector<double, VectFull>& correlation,
                    Vector<T1, Collection, Allocator1>& output);

    };


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_BASEPERTURBATIONMANAGER_HXX
#endif
