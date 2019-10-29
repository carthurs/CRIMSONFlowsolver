// Copyright (C) 2010, INRIA
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


#ifndef VERDANDI_FILE_SHARE_VERDANDIOPS_HXX


#include "ops/OpsHeader.hxx"


namespace Verdandi
{


    //! This class extends the Ops::Ops class with Verdandi-related features.
    class VerdandiOps: public ::Ops::Ops
    {
    public:
        // Constructor and destructor.
        VerdandiOps();
        VerdandiOps(string file_path);
        ~VerdandiOps();

    public:
        template<class TD, class T>
        void Set(string name, string constraint, const TD& default_value,
                 T& value);
        template<class T>
        void Set(string name, string constraint, T& value);
        template <class T>
        void Set(string name, T& value);

        template<class T>
        T Get(string name);
        template<class T>
        T Get(string name, string constraint);
        template<class T>
        T Get(string name, string constraint, const T& default_value);

        template<class T>
        bool Is(string name);

    protected:
        using ::Ops::Ops::SetValue;
        template<class T, class Allocator>
        void SetValue(string name, string constraint,
                      const Seldon::Vector<T, VectFull, Allocator>&
                      default_value,
                      bool with_default,
                      Seldon::Vector<T, VectFull, Allocator>& value);
        template<class T, class Allocator>
        void SetValue(string name, string constraint,
                      const vector<Seldon::Vector<T, VectFull, Allocator> >&
                      default_value,
                      bool with_default,
                      vector<Seldon::Vector<T, VectFull, Allocator> >& value);
        template<class T, class Prop, class Storage, class Allocator>
        void SetValue(string name, string constraint,
                      const Seldon::Matrix<T, Prop, Storage, Allocator>&
                      default_value,
                      bool with_default,
                      Seldon::Matrix<T, Prop, Storage, Allocator>& value);
        template<class T, class Prop, class Storage, class Allocator>
        void SetValue(string name, string constraint,
                      const
                      vector<Seldon::Matrix<T, Prop, Storage, Allocator> >&
                      default_value,
                      bool with_default,
                      vector<Seldon::Matrix<T, Prop, Storage, Allocator> >&
                      value);

        using ::Ops::Ops::IsParam;
        template<class T, class Allocator>
        bool IsParam(string name,
                     Seldon::Vector<T, VectFull, Allocator>& value);
        template<class T, class Prop, class Storage, class Allocator>
        bool IsParam(string name,
                     Seldon::Matrix<T, Prop, Storage, Allocator>& value);

        using ::Ops::Ops::CheckConstraint;
        bool CheckConstraint(string expression);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_SHARE_VERDANDIOPS_HXX
#endif
