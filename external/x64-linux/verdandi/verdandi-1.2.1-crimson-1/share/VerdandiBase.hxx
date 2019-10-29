// Copyright (C) 2009 INRIA
// Author(s): Claire Mouton
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


#ifndef VERDANDI_FILE_SHARE_VERDANDIBASE_HXX


namespace Verdandi
{


    //////////////////
    // VERDANDIBASE //
    //////////////////


    //! Base class for Verdandi objects.
    /*!
      \tparam T the type of floating-point numbers.
    */
    class VerdandiBase
    {
    public:
        // Constructor and destructor.
        VerdandiBase();
        virtual ~VerdandiBase();

        virtual string GetName() const;
        virtual void Message(string message);
        static void StaticMessage(void* object, string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_SHARE_VERDANDIBASE_HXX
#endif
