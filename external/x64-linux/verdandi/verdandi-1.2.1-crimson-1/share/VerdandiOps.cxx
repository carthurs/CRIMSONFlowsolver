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


#ifndef VERDANDI_FILE_SHARE_VERDANDIOPS_CXX

#include "ops/Ops.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    /*! Nothing is performed. The Lua state is set to NULL.
     */
    VerdandiOps::VerdandiOps(): ::Ops::Ops()
    {
    }


    //! Main constructor.
    /*! The Lua configuration file is loaded and run. An exception may be
      raised during this evaluation.
      \param[in] file_path path to the configuration file.
    */
    VerdandiOps::VerdandiOps(string file_path): ::Ops::Ops(file_path)
    {
    }


    //! Destructor.
    /*! Destroys the Lua state object.
     */
    VerdandiOps::~VerdandiOps()
    {
    }


    //////////////////
    // MAIN METHODS //
    //////////////////


    //! Retrieves a value from the configuration file.
    /*!
      \param[in] name name of the entry.
      \param[in] constraint constraint that the entry value must satisfy.
      \param[in] default_value default value for the entry in case it is not
      found in the configuration file.
      \param[out] value value of the entry.
    */
    template<class TD, class T>
    void
    VerdandiOps::Set(string name, string constraint,
                     const TD& default_value, T& value)
    {
        SetValue(name, constraint, default_value, true, value);
    }


    //! Retrieves a value from the configuration file.
    /*!
      \param[in] name name of the entry.
      \param[in] constraint constraint that the entry value must satisfy.
      \param[out] value value of the entry.
    */
    template<class T>
    void VerdandiOps::Set(string name, string constraint, T& value)
    {
        SetValue(name, constraint, value, false, value);
    }


    //! Retrieves a value from the configuration file.
    /*!
      \param[in] name name of the entry.
      \param[out] value value of the entry.
    */
    template <class T>
    void VerdandiOps::Set(string name, T& value)
    {
        SetValue(name, "", value, false, value);
    }


    //! Retrieves a value from the configuration file.
    /*!
      \param[in] name name of the entry.
      \param[in] constraint constraint that the entry value must satisfy.
      \param[in] default_value default value for the entry in case it is not
      found in the configuration file.
      \return The value of the entry.
    */
    template<class T>
    T VerdandiOps::Get(string name, string constraint, const T& default_value)
    {
        T value;
        SetValue(name, constraint, default_value, true, value);
        return value;
    }


    //! Retrieves a value from the configuration file.
    /*!
      \param[in] name name of the entry.
      \param[in] constraint constraint that the entry value must satisfy.
      \return The value of the entry.
    */
    template<class T>
    T VerdandiOps::Get(string name, string constraint)
    {
        T value;
        SetValue(name, constraint, value, false, value);
        return value;
    }


    //! Retrieves a value from the configuration file.
    /*!
      \param[in] name name of the entry.
      \return The value of the entry.
    */
    template <class T>
    T VerdandiOps::Get(string name)
    {
        T value;
        SetValue(name, "", value, false, value);
        return value;
    }


    //! Checks whether \a name is of type 'T'.
    /*! On exit, the value of the entry (if it exists) is on the stack.
      \param[in] name the name of the entry whose type is checked.
      \return True if the entry is of type 'T', false otherwise.
      \note The prefix is prepended to \a name. If \a name does not exist, an
      exception is raised.
    */
    template<class T>
    bool VerdandiOps::Is(string name)
    {
        T value;
        return IsParam(name, value);
    }


    //! Retrieves a value and checks if it satisfies given constraints.
    /*! If the entry is not found, the default value is returned (if any). If
      the vector is given in a file, the Lua entry should provide the file
      name. The file is read with method 'Read' of \a value.
      \param[in] name name of the entry.
      \param[in] constraint constraint to be satisfied.
      \param[in] default_value default value.
      \param[in] with_default is there a default value? If not, \a
      default_value is ignored.
      \param[out] value the value of the entry named \a name.
      \note The default value may not satisfy the constraint.
    */
    template<class T, class Allocator>
    void VerdandiOps::SetValue(string name, string constraint,
                               const Seldon::Vector<T, VectFull, Allocator>&
                               default_value,
                               bool with_default,
                               Seldon::Vector<T, VectFull, Allocator>& value)
    {
        if (!this->Exists(name))
        {
            if (with_default)
                value = default_value;
            else
                throw Error("SetValue",
                            "The " + Entry(name) + " was not found.");
        }
        else if (this->Is<string>(name))
        {
            string filename;
            SetValue(name, "", "", false, filename);
            value.Read(filename);
            if (constraint.empty())
                return;
            for (int i = 0; i < value.GetLength(); i++)
                if (!CheckConstraintOnValue(to_str(value(i)), constraint))
                    throw Error("SetValue",
                                "The entry "
                                + Entry(name + "[" + to_str(i + 1) + "]")
                                + " does not satisfy the constraint:\n"
                                + Constraint(constraint));
        }
        else
        {
            std::vector<T> default_data(default_value.GetLength()), data;
            for (int i = 0; i < default_value.GetLength(); i++)
                default_data[i] = default_value(i);
            SetValue(name, constraint, default_data, with_default, data);
            value.Reallocate(int(data.size()));
            for (size_t i = 0; i < data.size(); i++)
                value(i) = data[i];
        }
    }


    //! Retrieves a value and checks if it satisfies given constraints.
    /*! If the entry is not found, the default value is returned (if any).

      If the list of vectors is given by a Lua array, it is assumed that the
      number of vectors is given by the length of \a value. If \a value is
      empty, it is assumed there is only one vector. All vectors in \a value
      should have the same length, otherwise an exception is raised.

      If the list of vectors is given in a file, the Lua entry should provide
      the file name. The file is read with method 'Read' of \a value.
      \param[in] name name of the entry.
      \param[in] constraint constraint to be satisfied.
      \param[in] default_value default value.
      \param[in] with_default is there a default value? If not, \a
      default_value is ignored.
      \param[out] value the value of the entry named \a name.
      \note The default value may not satisfy the constraint.
    */
    template<class T, class Allocator>
    void VerdandiOps
    ::SetValue(string name, string constraint,
               const vector<Seldon::Vector<T, VectFull, Allocator> >&
               default_value,
               bool with_default,
               vector<Seldon::Vector<T, VectFull, Allocator> >& value)
    {
        if (!this->Exists(name))
        {
            if (with_default)
                value = default_value;
            else
                throw Error("SetValue",
                            "The " + Entry(name) + " was not found.");
        }
        else if (this->Is<string>(name))
        {
            string filename;
            value.clear();
            SetValue(name, "", "", false, filename);
            ifstream f(filename.c_str());
            while (f.peek() != -1)
            {
                value.push_back(Seldon::Vector<T, VectFull, Allocator>());
                value[value.size() - 1].Read(f);
            }
            if (constraint.empty())
                return;
            for (size_t i = 0; i < value.size(); i++)
                for (int j = 0; j < value[i].GetLength(); j++)
                    if (!CheckConstraintOnValue(to_str(value[i](j)),
                                                constraint))
                        throw Error("SetValue",
                                    "The entry "
                                    + Entry(name + "[" + to_str(i + 1) + "]"
                                            + "[" + to_str(j + 1) + "]")
                                    + " does not satisfy the constraint:\n"
                                    + Constraint(constraint));
        }
        else
        {
            int N = 1;
            for (size_t i = 0; i < default_value.size(); i++)
                N += default_value[i].GetLength();
            std::vector<T> default_data(N);
            int index = 0;
            for (size_t i = 0; i < default_value.size(); i++)
                for (int j = 0; j < default_value[i].GetLength(); j++)
                    default_data[index++] = default_value[i](j);

            std::vector<T> data;
            SetValue(name, constraint, default_data, with_default, data);
            if (value.size() != 0 && data.size() % value.size() != 0)
                throw Error("SetValue",
                            "The entry " + Entry(name)
                            + " contains " + to_str(data.size())
                            + " elements, which is incompatible with a list"
                            + " of " + to_str(value.size()) + " vectors of "
                            "the same length.");
            if (value.size() == 0)
                value.resize(1);
            index = 0;
            for (size_t i = 0; i < value.size(); i++)
                value[i].Reallocate(int(data.size() / value.size()));
            for (size_t i = 0; i < value.size(); i++)
                for (int j = 0; j < value[i].GetLength(); j++)
                    value[i](j) = data[index++];
        }
    }


    //! Retrieves a value and checks if it satisfies given constraints.
    /*! If the entry is not found, the default value is returned (if any).

      If the matrix is given by a Lua array, it is assumed that the Lua array
      is in row-major order (the rows are stored one after the other), and the
      number of rows is given by the matrix on entry. If the matrix on entry
      has no row, it is assumed the matrix to be read has a single row. The
      matrix is reallocated according to the number of columns found in the
      Lua array. If the number of elements in the Lua array is not a multiple
      of the number of rows in the matrix, an exception is raised.

      If the matrix is given in a file, the Lua entry should provide the file
      name. The file is read with method 'Read' of \a value.
      \param[in] name name of the entry.
      \param[in] constraint constraint to be satisfied.
      \param[in] default_value default value.
      \param[in] with_default is there a default value? If not, \a
      default_value is ignored.
      \param[out] value the value of the entry named \a name.
      \note The default value may not satisfy the constraint.
    */
    template<class T, class Prop, class Storage, class Allocator>
    void VerdandiOps
    ::SetValue(string name, string constraint,
               const Seldon::Matrix<T, Prop, Storage, Allocator>&
               default_value,
               bool with_default,
               Seldon::Matrix<T, Prop, Storage, Allocator>& value)
    {
        if (!this->Exists(name))
        {
            if (with_default)
                value = default_value;
            else
                throw Error("SetValue",
                            "The " + Entry(name) + " was not found.");
        }
        else if (this->Is<string>(name))
        {
            string filename;
            SetValue(name, "", "", false, filename);
            value.Read(filename);
            if (constraint.empty())
                return;
            for (int i = 0; i < value.GetM(); i++)
                for (int j = 0; j < value.GetN(); j++)
                    if (!CheckConstraintOnValue(to_str(value(i, j)),
                                                constraint))
                        throw Error("SetValue",
                                    "The entry "
                                    + Entry(name + "[" + to_str(i + 1) + ", "
                                            + to_str(j + 1) + "]")
                                    + " does not satisfy the constraint:\n"
                                    + Constraint(constraint));
        }
        else
        {

            std::vector<T> default_data(default_value.GetM()
                                        * default_value.GetN());
            int index = 0;
            for (int i = 0; i < default_value.GetM(); i++)
                for (int j = 0; j < default_value.GetN(); j++)
                    default_data[index++] = default_value(i, j);

            std::vector<T> data;
            SetValue(name, constraint, default_data, with_default, data);
            if (value.GetM() != 0 && int(data.size()) % value.GetM() != 0)
                throw Error("SetValue",
                            "The entry " + Entry(name)
                            + " contains " + to_str(int(data.size()))
                            + " elements, which is incompatible with a matrix"
                            + " containing " + to_str(value.GetM())
                            + " rows.");
            if (value.GetM() == 0)
                value.Reallocate(1, int(data.size()));
            else
                value.Reallocate(value.GetM(),
                                 int(data.size()) / value.GetM());
            index = 0;
            for (int i = 0; i < value.GetM(); i++)
                for (int j = 0; j < value.GetN(); j++)
                    value(i, j) = data[index++];
        }
    }


    //! Retrieves a value and checks if it satisfies given constraints.
    /*! If the entry is not found, the default value is returned (if any).

      If the list of matrices is given by a Lua array, it is assumed that (1)
      the number of matrices is given by the length of \a value; (2) the
      number of rows in the matrices is given by the number of rows in the
      first matrix of \a value; (3) the Lua array is in row-major order (the
      rows are stored one after the other). If \a value is empty, it is filled
      with a single matrix with one row. If \a value is not empty but its
      first matrix is, the matrices are supposed to have a single row on
      exit. The matrices are finally reallocated according to the number of
      columns found in the Lua array. If the number of elements in the Lua
      array is not a multiple of the number matrices times the number of rows
      in the matrices, an exception is raised.

      If the list of matrices is given in a file, the Lua entry should provide
      the file name. The file is read with method 'Read' of the matrices.
      \param[in] name name of the entry.
      \param[in] constraint constraint to be satisfied.
      \param[in] default_value default value.
      \param[in] with_default is there a default value? If not, \a
      default_value is ignored.
      \param[out] value the value of the entry named \a name.
      \note The default value may not satisfy the constraint.
    */
    template<class T, class Prop, class Storage, class Allocator>
    void VerdandiOps
    ::SetValue(string name, string constraint,
               const
               vector<Seldon::Matrix<T, Prop, Storage, Allocator> >&
               default_value,
               bool with_default,
               vector<Seldon::Matrix<T, Prop, Storage, Allocator> >&
               value)
    {
        if (!this->Exists(name))
        {
            if (with_default)
                value = default_value;
            else
                throw Error("SetValue",
                            "The " + Entry(name) + " was not found.");
        }
        else if (this->Is<string>(name))
        {
            string filename;
            value.clear();
            SetValue(name, "", "", false, filename);
            ifstream f(filename.c_str());
            while (f.peek() != -1)
            {
                value
                    .push_back(Seldon::Matrix<T, Prop, Storage, Allocator>());
                value[value.size() - 1].Read(f);
            }
            if (constraint.empty())
                return;
            for (size_t i = 0; i < value.size(); i++)
                for (int j = 0; j < value[i].GetM(); j++)
                    for (int k = 0; k < value[i].GetN(); k++)
                        if (!CheckConstraintOnValue(to_str(value[i](j, k)),
                                                    constraint))
                            throw Error("SetValue",
                                        "The entry "
                                        + Entry(name + "[" + to_str(i + 1)
                                                + "]" + "[" + to_str(j + 1)
                                                + ", " + to_str(k + 1) + "]")
                                        + " does not satisfy the constraint:"
                                        + "\n" + Constraint(constraint));
        }
        else
        {
            int N = 1;
            for (size_t i = 0; i < default_value.size(); i++)
                N += default_value[i].GetM() * default_value[i].GetN();
            std::vector<T> default_data(N);
            int index = 0;
            for (size_t i = 0; i < default_value.size(); i++)
                for (int j = 0; j < default_value[i].GetM(); j++)
                    for (int k = 0; k < default_value[i].GetN(); k++)
                        default_data[index++] = default_value[i](j, k);

            std::vector<T> data;
            SetValue(name, constraint, default_data, with_default, data);

            int m = 1;
            if (value.size() != 0 && value[0].GetM() != 0)
                m = value[0].GetM();
            if (value.size() != 0 && data.size() % (m * value.size()) != 0)
                throw Error("SetValue",
                            "The entry " + Entry(name)
                            + " contains " + to_str(data.size())
                            + " elements, which is incompatible with a list"
                            + " of " + to_str(value.size()) + " matrices"
                            + " with " + to_str(m) + " rows.");
            if (value.size() == 0)
                value.resize(1);
            index = 0;
            int n = int(data.size() / value.size()) / m;
            for (size_t i = 0; i < value.size(); i++)
                value[i].Reallocate(m, n);
            for (size_t i = 0; i < value.size(); i++)
                for (int j = 0; j < value[i].GetM(); j++)
                    for (int k = 0; k < value[i].GetN(); k++)
                        value[i](j, k) = data[index++];
        }
    }


    //! Checks whether \a name is a 'Vector<T>'.
    /*!
      \param[in] name the name of the entry whose type is checked.
      \param[in] value anything: it is used to determine the type.
      \return True if the entry is a 'Vector<T>', false otherwise.
      \note The prefix is prepended to \a name. If \a name does not exist, an
      exception is raised.
    */
    template<class T, class Allocator>
    bool VerdandiOps::IsParam(string name,
                              Seldon::Vector<T, VectFull, Allocator>& value)
    {
        string str_value;
        std::vector<T> vect_value;
        return this->IsParam(name, str_value)
            || this->IsParam(name, vect_value);
    }


    //! Checks whether \a name is a 'Matrix<T>'.
    /*!
      \param[in] name the name of the entry whose type is checked.
      \param[in] value anything: it is used to determine the type.
      \return True if the entry is a 'Matrix<T>', false otherwise.
      \note The prefix is prepended to \a name. If \a name does not exist, an
      exception is raised.
    */
    template<class T, class Prop, class Storage, class Allocator>
    bool VerdandiOps
    ::IsParam(string name,
              Seldon::Matrix<T, Prop, Storage, Allocator>& value)
    {
        string str_value;
        std::vector<T> vect_value;
        return this->IsParam(name, str_value)
            || this->IsParam(name, vect_value);
    }


} // namespace Verdandi.


#define VERDANDI_FILE_SHARE_VERDANDIOPS_CXX
#endif

