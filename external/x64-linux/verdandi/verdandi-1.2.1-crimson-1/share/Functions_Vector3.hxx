// Copyright (C) 2010-2012 INRIA
// Author(s): Vivien Mallet, Anne Tilloy
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


#ifndef VERDANDI_FILE_SHARE_FUNCTIONS_VECTOR3_HXX


namespace Verdandi
{


    template <class T,
              class TV, class Allocator0, class Allocator1, class Allocator2>
    void RemoveData(T value, Vector3<TV, Allocator0, Allocator1,
                    Allocator2>& V);


} // namespace Verdandi


#define VERDANDI_FILE_SHARE_FUNCTIONS_VECTOR3_HXX
#endif
