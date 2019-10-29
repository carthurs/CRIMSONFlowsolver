// Copyright (C) 2010, Vivien Mallet
//
// This file is part of Ops, a library for parsing Lua configuration files.
//
// Ops is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Ops is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Ops. If not, see http://www.gnu.org/licenses/.


#ifndef OPS_FILE_CLASSOPS_CXX


#include "ClassOps.hxx"


namespace Ops
{


  /////////////////////////////////
  // CONSTRUCTORS AND DESTRUCTOR //
  /////////////////////////////////


  //! Default constructor.
  /*! Nothing is performed. A Lua state is opened.
   */
  Ops::Ops()
  {
    state_ = lua_open();
    luaL_openlibs(state_);

    // Defines 'ops_in' for the user. It checks whether an element is in a
    // table.
    string code = "function ops_in(v, table)\
    for _, value in ipairs(table) do        \
        if v == value then                  \
            return true                     \
        end                                 \
    end                                     \
    return false                            \
    end";
    if (luaL_dostring(state_, code.c_str()))
      throw Error("Ops()", lua_tostring(state_, -1));
  }


  //! Main constructor.
  /*! The Lua configuration file is loaded and run. An exception may be raised
    during this evaluation.
    \param[in] file_path path to the configuration file.
  */
  Ops::Ops(string file_path):
    file_path_(file_path), state_(NULL)
  {
    Open(file_path_);
  }


  //! Destructor.
  /*! Destroys the Lua state object.
   */
  Ops::~Ops()
  {
    Close();
  }


  //////////////////
  // MAIN METHODS //
  //////////////////


  //! Opens a new configuration file.
  /*! The previous configuration file (if any) is closed. The prefix is
    cleared.
    \param[in] file_path path to the configuration file.
    \param[in] close_state should the Lua state be closed?
  */
  void Ops::Open(string file_path, bool close_state)
  {
    if (close_state)
      {
        Close();
        state_ = lua_open();
        luaL_openlibs(state_);
        // Defines 'ops_in' for the user. It checks whether an element is in a
        // table.
        string code = "function ops_in(v, table)\
        for _, value in ipairs(table) do        \
            if v == value then                  \
                return true                     \
            end                                 \
        end                                     \
        return false                            \
        end";
        if (luaL_dostring(state_, code.c_str()))
          throw Error("Open(string, bool)", lua_tostring(state_, -1));
      }

    ClearPrefix();
    file_path_ = file_path;
    if (luaL_dofile(state_, file_path_.c_str()))
      throw Error("Open(string, bool)", lua_tostring(state_, -1));
  }


  //! Reloads the current configuration file.
  /*! The configuration file is closed and reopened.
    \param[in] close_state should the Lua state be closed before reloading the
    file?
  */
  void Ops::Reload(bool close_state)
  {
    Open(file_path_, close_state);
  }


  //! Closes the configuration file (if any is open).
  /*! Destroys the Lua state object. The prefix is cleared.
    \param[in] file_path path to the configuration file.
  */
  void Ops::Close()
  {
    ClearPrefix();
    read_bool.clear();
    read_int.clear();
    read_float.clear();
    read_double.clear();
    read_string.clear();
    read_vect_bool.clear();
    read_vect_int.clear();
    read_vect_float.clear();
    read_vect_double.clear();
    read_vect_string.clear();
    if (state_ != NULL)
      lua_close(state_);
    state_ = NULL;
  }


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
  Ops::Set(string name, string constraint, const TD& default_value, T& value)
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
  void Ops::Set(string name, string constraint, T& value)
  {
    SetValue(name, constraint, value, false, value);
  }


  //! Retrieves a value from the configuration file.
  /*!
    \param[in] name name of the entry.
    \param[out] value value of the entry.
  */
  template <class T>
  void Ops::Set(string name, T& value)
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
  T Ops::Get(string name, string constraint, const T& default_value)
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
  T Ops::Get(string name, string constraint)
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
  T Ops::Get(string name)
  {
    T value;
    SetValue(name, "", value, false, value);
    return value;
  }


  //! Applies a Lua function.
  /*!
    \param[in] name name of the function.
    \param[in] in parameters of the function.
    \param[out] out outputs of the function.
    \note The prefix is prepended to \a name.
  */
  template<class Tin, class Tout>
  void Ops::Apply(string name, const std::vector<Tin>& in,
                  std::vector<Tout>& out)
  {
    PutOnStack(Name(name));
    PushOnStack(in);

    int n = lua_gettop(state_);

    if (lua_pcall(state_, int(in.size()), LUA_MULTRET, 0) != 0)
      throw Error("Apply(string, vector, vector&)",
                  "While calling " + Function(name) + ":\n  "
                  + lua_tostring(state_, -1));

    n = lua_gettop(state_) - n + 1 + int(in.size());

    out.resize(n);
    for (int i = 0; i < n; i++)
      if (!Convert(i - n, out[i]))
        {
          std::ostringstream str;
          str << i;
          throw Error("Apply(string, vector, vector&)",
                      "The returned value #" + str.str() + " of \""
                      + Name(name) + "\" is not of correct type.");
        }

    ClearStack();
  }


  //! Applies a Lua function.
  /*!
    \param[in] name name of the function.
    \param[in] arg0 parameter of the function.
    \return First output of the function.
    \note The prefix is prepended to \a name.
  */
  template<class T>
  T Ops::Apply(string name, const T& arg0)
  {
    std::vector<T> in, out;
    in.push_back(arg0);
    Apply(name, in, out);
    return out[0];
  }


  //! Applies a Lua function.
  /*!
    \param[in] name name of the function.
    \param[in] arg0 parameter of the function.
    \param[in] arg1 parameter of the function.
    \return First output of the function.
    \note The prefix is prepended to \a name.
  */
  template<class T>
  T Ops::Apply(string name, const T& arg0, const T& arg1)
  {
    std::vector<T> in, out;
    in.push_back(arg0);
    in.push_back(arg1);
    Apply(name, in, out);
    return out[0];
  }


  //! Applies a Lua function.
  /*!
    \param[in] name name of the function.
    \param[in] arg0 parameter of the function.
    \param[in] arg1 parameter of the function.
    \param[in] arg2 parameter of the function.
    \return First output of the function.
    \note The prefix is prepended to \a name.
  */
  template<class T>
  T Ops::Apply(string name, const T& arg0, const T& arg1, const T& arg2)
  {
    std::vector<T> in, out;
    in.push_back(arg0);
    in.push_back(arg1);
    in.push_back(arg2);
    Apply(name, in, out);
    return out[0];
  }


  //! Applies a Lua function.
  /*!
    \param[in] name name of the function.
    \param[in] arg0 parameter of the function.
    \param[in] arg1 parameter of the function.
    \param[in] arg2 parameter of the function.
    \param[in] arg3 parameter of the function.
    \return First output of the function.
    \note The prefix is prepended to \a name.
  */
  template<class T>
  T Ops::Apply(string name, const T& arg0, const T& arg1, const T& arg2,
               const T& arg3)
  {
    std::vector<T> in, out;
    in.push_back(arg0);
    in.push_back(arg1);
    in.push_back(arg2);
    in.push_back(arg3);
    Apply(name, in, out);
    return out[0];
  }


  //! Applies a Lua function.
  /*!
    \param[in] name name of the function.
    \param[in] arg0 parameter of the function.
    \param[in] arg1 parameter of the function.
    \param[in] arg2 parameter of the function.
    \param[in] arg3 parameter of the function.
    \param[in] arg4 parameter of the function.
    \return First output of the function.
    \note The prefix is prepended to \a name.
  */
  template<class T>
  T Ops::Apply(string name, const T& arg0, const T& arg1, const T& arg2,
               const T& arg3, const T& arg4)
  {
    std::vector<T> in, out;
    in.push_back(arg0);
    in.push_back(arg1);
    in.push_back(arg2);
    in.push_back(arg3);
    in.push_back(arg4);
    Apply(name, in, out);
    return out[0];
  }


  //! Returns the list of entries inside an entry.
  /*!
    \param[in] name name of the entry to search in.
    \return The list of entries under \a name, sorted in alphabetical order
    (with numbers coming first).
    \note The prefix is prepended to \a name before the search. If the entry
    \a name does not exist or does not contain other entries, an exception is
    raised.
  */
  std::vector<string> Ops::GetEntryList(string name)
  {
    PutOnStack(Name(name));

    if (lua_isnil(state_, -1))
      throw Error("GetEntryList",
                  "The " + Entry(name) + " was not found.");

    if (!lua_istable(state_, -1))
      throw Error("GetEntryList",
                  "The " + Entry(name) + " does not contain other entries.");

    std::vector<string> key_list;
    string key;
    // Now loops over all elements of the table.
    lua_pushnil(state_);
    while (lua_next(state_, -2) != 0)
      {
        // Duplicates the key so that 'lua_tostring' (applied to it) should
        // not interfere with 'lua_next'.
        lua_pushvalue(state_, -2);

        if (!Convert(-1, key))
          throw Error("GetEntryList",
                      "Unable to read the keys of " + Entry(name) + ".");
        key_list.push_back(key);

        lua_pop(state_, 2);
      }

    // Sorts the keys.
    std::sort(key_list.begin(), key_list.end());

    ClearStack();

    return key_list;
  }


  //! Checks that a certain entry satisfies a constraint.
  /*!
    \param[in] name the name of the entry whose consistency with \a constraint
    is to be checked.
    \param[in] constraint the constraint to be satisfied.
    \return True if the constraint is satisfied, false otherwise.
  */
  bool Ops::CheckConstraint(string name, string constraint)
  {
    if (constraint == "")
      return true;

    string code;
    code = "function ops_check_constraint(v)\nreturn " + constraint \
      + "\nend\nops_result = ops_check_constraint(" + Name(name) + ")";
    if (luaL_dostring(state_, code.c_str()))
      throw Error("CheckConstraint",
                  "While checking " + Entry(name) + ":\n  "
                  + string(lua_tostring(state_, -1)));

    PutOnStack("ops_result");
    if (!lua_isboolean(state_, -1))
      throw Error("CheckConstraint",
                  "For " + Entry(name) + ", the following constraint did "
                  "not return a Boolean:\n" + Constraint(constraint));

    return static_cast<bool>(lua_toboolean(state_, -1));
  }


  //! Checks that a certain value satisfies a constraint.
  /*!
    \param[in] value value whose consistency with \a constraint is to be
    checked.
    \param[in] constraint the constraint to be satisfied.
    \return True if the constraint is satisfied, false otherwise.
  */
  bool Ops::CheckConstraintOnValue(string value, string constraint)
  {
    if (constraint == "")
      return true;

    string code;
    code = "function ops_check_constraint(v)\nreturn " + constraint \
      + "\nend\nops_result = ops_check_constraint(" + value + ")";
    if (luaL_dostring(state_, code.c_str()))
      throw Error("CheckConstraintOnValue",
                  "While checking the value \"" + value + "\":\n  "
                  + string(lua_tostring(state_, -1)));

    PutOnStack("ops_result");
    if (!lua_isboolean(state_, -1))
      throw Error("CheckConstraint",
                  "For value \"" + value + "\", the following constraint did "
                  "not return a Boolean:\n" + Constraint(constraint));

    return static_cast<bool>(lua_toboolean(state_, -1));
  }


  //! Puts \a name on top of the stack.
  /*! If \a name is a simple variable, it calls 'lua_getglobal' once. But if
    \a name is encapsulated in a table, this method iterates until it finds
    the variable.
    \param[in] name the name of the entry to be put on top of the stack.
    \note The prefix is not prepended to \a name.
  */
  void Ops::PutOnStack(string name)
  {
    if (name.empty())
      {
        lua_pushvalue(state_, LUA_GLOBALSINDEX);
        return;
      }

    size_t end = name.find_first_of(".[");
    if (end == 0)
      // The name starts with '.' or '[': wrong syntax.
      {
        lua_pushnil(state_);
        return;
      }
    if (end == string::npos)
      {
        lua_getglobal(state_, name.c_str());
        return;
      }

    lua_getglobal(state_, name.substr(0, end).c_str());

    if (name[end] == '.')
      WalkDown(name.substr(end + 1).c_str());
    else
      WalkDown(name.substr(end).c_str());
  }


  //! Checks whether \a name exists.
  /*! On exit, the value of the entry (if it exists) is on the stack.
    \param[in] name the name of the entry whose existence is checked.
    \return True if the entry exists, false otherwise.
    \note The prefix is prepended to \a name.
  */
  bool Ops::Exists(string name)
  {
    PutOnStack(Name(name));
    bool exists = !lua_isnil(state_, -1);
    ClearStack();
    return exists;
  }


  //! Checks whether \a name is of type 'T'.
  /*! On exit, the value of the entry (if it exists) is on the stack.
    \param[in] name the name of the entry whose type is checked.
    \return True if the entry is of type 'T', false otherwise.
    \note The prefix is prepended to \a name. If \a name does not exist, an
    exception is raised.
  */
  template<class T>
  bool Ops::Is(string name)
  {
    T value;
    return IsParam(name, value);
  }


  //! Checks whether \a name is a table.
  /*! On exit, the value of the entry (if it exists) is on the stack.
    \param[in] name the name of the entry whose type is checked.
    \return True if the entry is a table, false otherwise.
    \note The prefix is prepended to \a name. If \a name does not exist, an
    exception is raised.
  */
  bool Ops::IsTable(string name)
  {
    PutOnStack(Name(name));
    return lua_istable(state_, -1);
  }


  //! Checks whether \a name is a function.
  /*! On exit, the value of the entry (if it exists) is on the stack.
    \param[in] name the name of the entry whose type is checked.
    \return True if the entry is a function, false otherwise.
    \note The prefix is prepended to \a name. If \a name does not exist, an
    exception is raised.
  */
  bool Ops::IsFunction(string name)
  {
    PutOnStack(Name(name));
    return lua_isfunction(state_, -1);
  }


  //! Pushes an element onto the stack.
  /*!
    \param[in] value element to be pushed.
  */
  void Ops::PushOnStack(bool value)
  {
    lua_pushboolean(state_, value);
  }


  //! Pushes an element onto the stack.
  /*!
    \param[in] value element to be pushed.
  */
  void Ops::PushOnStack(int value)
  {
    lua_pushinteger(state_, value);
  }


  //! Pushes an element onto the stack.
  /*!
    \param[in] value element to be pushed.
  */
  void Ops::PushOnStack(float value)
  {
    lua_pushnumber(state_, value);
  }


  //! Pushes an element onto the stack.
  /*!
    \param[in] value element to be pushed.
  */
  void Ops::PushOnStack(double value)
  {
    lua_pushnumber(state_, value);
  }


  //! Pushes an element onto the stack.
  /*!
    \param[in] value element to be pushed.
  */
  void Ops::PushOnStack(string value)
  {
    lua_pushstring(state_, value.c_str());
  }


  //! Pushes a vector onto the stack.
  /*! Every element of the vector \a v is pushed onto the stack, in the same
    order as in the vector.
    \param[in] v vector to be pushed.
  */
  template<class T>
  void Ops::PushOnStack(const std::vector<T>& v)
  {
    for (int i = 0; i < int(v.size()); i++)
      PushOnStack(v[i]);
  }


  //! Clears the stack.
  void Ops::ClearStack()
  {
    lua_pop(state_, lua_gettop(state_));
  }


  //! Execute Lua code.
  /*!
    \param[in] file_path path to the file to be processed.
  */
  void Ops::DoFile(string file_path)
  {
    if (luaL_dofile(state_, file_path.c_str()))
      throw Error("DoFile(string)", lua_tostring(state_, -1));
  }


  //! Execute Lua code.
  /*!
    \param[in] expression the Lua code to be evaluated.
  */
  void Ops::DoString(string expression)
  {
    if (luaL_dostring(state_, expression.c_str()))
      throw Error("DoString(string)", lua_tostring(state_, -1));
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////


  //! Returns the path to the configuration file.
  /*!
    \return The path to the configuration file.
  */
  string Ops::GetFilePath() const
  {
    return file_path_;
  }


  //! Returns the Lua state object.
  /*!
    \return The Lua state object.
  */
  lua_State* Ops::GetState()
  {
    return state_;
  }


  //! Returns the Lua state object.
  /*!
    \return The Lua state object.
  */
  const lua_State* Ops::GetState() const
  {
    return state_;
  }


  //! Returns the current prefix.
  /*!
    \return The current prefix.
  */
  string Ops::GetPrefix() const
  {
    return prefix_;
  }


  //! Sets the prefix.
  /*!
    \param[in] prefix the new prefix.
  */
  void Ops::SetPrefix(string prefix)
  {
    prefix_ = prefix;
  }


  //! Clears the current prefix.
  void Ops::ClearPrefix()
  {
    prefix_ = "";
  }


  //! Returns the list of the entry names that were read.
  /*!
    \return The list of the entries that were read, except the functions. The
    names are sorted.
  */
  std::vector<string> Ops::GetReadEntryList()
  {
    std::vector<string> name_list;
    AppendKey(read_bool, name_list);
    AppendKey(read_int, name_list);
    AppendKey(read_float, name_list);
    AppendKey(read_double, name_list);
    AppendKey(read_string, name_list);
    AppendKey(read_vect_bool, name_list);
    AppendKey(read_vect_int, name_list);
    AppendKey(read_vect_float, name_list);
    AppendKey(read_vect_double, name_list);
    AppendKey(read_vect_string, name_list);

    sort(name_list.begin(), name_list.end());
    name_list.erase(std::unique(name_list.begin(), name_list.end()),
                    name_list.end());

    return name_list;
  }


  //! Updates the values of all read variables.
  /*! After a variable is read, it can be modified with calls to 'DoFile' or
    'DoString'. For 'LuaDefinition' and 'WriteLuaDefinition' to use the
    updated value, this method should be called. It will read all variables
    again.
  */
  void Ops::UpdateLuaDefinition()
  {
    string prefix = prefix_;

    for (std::map<string, bool>::iterator i = read_bool.begin();
         i != read_bool.end(); i++)
      Get(i->first, "", i->second);
    for (std::map<string, int>::iterator i = read_int.begin();
         i != read_int.end(); i++)
      Get(i->first, "", i->second);
    for (std::map<string, float>::iterator i = read_float.begin();
         i != read_float.end(); i++)
      Get(i->first, "", i->second);
    for (std::map<string, double>::iterator i = read_double.begin();
         i != read_double.end(); i++)
      Get(i->first, "", i->second);
    for (std::map<string, string>::iterator i = read_string.begin();
         i != read_string.end(); i++)
      Get(i->first, "", i->second);

    for (std::map<string, std::vector<bool> >::iterator
           i = read_vect_bool.begin();
         i != read_vect_bool.end(); i++)
      Get(i->first, "", i->second);
    for (std::map<string, std::vector<int> >::iterator
           i = read_vect_int.begin();
         i != read_vect_int.end(); i++)
      Get(i->first, "", i->second);
    for (std::map<string, std::vector<float> >::iterator
           i = read_vect_float.begin();
         i != read_vect_float.end(); i++)
      Get(i->first, "", i->second);
    for (std::map<string, std::vector<double> >::iterator
           i = read_vect_double.begin();
         i != read_vect_double.end(); i++)
      Get(i->first, "", i->second);
    for (std::map<string, std::vector<string> >::iterator
           i = read_vect_string.begin();
         i != read_vect_string.end(); i++)
      Get(i->first, "", i->second);

    prefix_ = prefix;
  }


  //! Returns a Lua line that defines an entry already read.
  /*! If the entry was not already read, an exception is thrown.
    \param[in] name name of the entry.
    \return A Lua line that defines \a name. It is in the form "name = value".
    \warning The entry may have any type supported by Ops, except a function.
    \note The prefix is not prepended to \a name.
  */
  string Ops::LuaDefinition(string name)
  {
    std::ostringstream output;

    // Below, the name is searched in all maps. When the name is found, its
    // value is converted to a string (that Lua can process).

    std::map<string, bool>::iterator i_bool = read_bool.find(name);
    if (i_bool != read_bool.end())
      if (i_bool->second)
        output << name << " = true";
      else
        output << name << " = false";
    // In case the entry was read under two different types, only one type is
    // returned. So, once the entry is found, this method returns the
    // definition.
    if (!output.str().empty())
      return output.str();

    std::map<string, int>::iterator i_int = read_int.find(name);
    if (i_int != read_int.end())
      output << name << " = " << i_int->second;
    if (!output.str().empty())
      return output.str();

    std::map<string, float>::iterator i_float = read_float.find(name);
    if (i_float != read_float.end())
      output << name << " = " << i_float->second;
    if (!output.str().empty())
      return output.str();

    std::map<string, double>::iterator i_double = read_double.find(name);
    if (i_double != read_double.end())
      output << name << " = " << i_double->second;
    if (!output.str().empty())
      return output.str();

    std::map<string, string>::iterator i_string = read_string.find(name);
    if (i_string != read_string.end())
      output << name << " = \"" << i_string->second << "\"";
    if (!output.str().empty())
      return output.str();

    std::map<string, std::vector<bool> >::iterator i_vect_bool
      = read_vect_bool.find(name);
    if (i_vect_bool != read_vect_bool.end())
      {
        output << name << " = {";
        int size = i_vect_bool->second.size();
        for (int i = 0; i < size - 1; i++)
          if (i_vect_bool->second[i])
            output << "true, ";
          else
            output << "false, ";
        if (size != 0)
          if (i_vect_bool->second[size - 1])
            output << "true";
          else
            output << "false";
        output << "}";
      }
    if (!output.str().empty())
      return output.str();

    std::map<string, std::vector<int> >::iterator i_vect_int
      = read_vect_int.find(name);
    if (i_vect_int != read_vect_int.end())
      {
        output << name << " = {";
        int size = i_vect_int->second.size();
        for (int i = 0; i < size - 1; i++)
          output << i_vect_int->second[i] << ", ";
        if (size != 0)
          output << i_vect_int->second[size - 1];
        output << "}";
      }
    if (!output.str().empty())
      return output.str();

    std::map<string, std::vector<float> >::iterator i_vect_float
      = read_vect_float.find(name);
    if (i_vect_float != read_vect_float.end())
      {
        output << name << " = {";
        int size = i_vect_float->second.size();
        for (int i = 0; i < size - 1; i++)
          output << i_vect_float->second[i] << ", ";
        if (size != 0)
          output << i_vect_float->second[size - 1];
        output << "}";
      }
    if (!output.str().empty())
      return output.str();

    std::map<string, std::vector<double> >::iterator i_vect_double
      = read_vect_double.find(name);
    if (i_vect_double != read_vect_double.end())
      {
        output << name << " = {";
        int size = i_vect_double->second.size();
        for (int i = 0; i < size - 1; i++)
          output << i_vect_double->second[i] << ", ";
        if (size != 0)
          output << i_vect_double->second[size - 1];
        output << "}";
      }
    if (!output.str().empty())
      return output.str();

    std::map<string, std::vector<string> >::iterator i_vect_string
      = read_vect_string.find(name);
    if (i_vect_string != read_vect_string.end())
      {
        output << name << " = {";
        int size = i_vect_string->second.size();
        for (int i = 0; i < size - 1; i++)
          output << "\"" << i_vect_string->second[i] << "\", ";
        if (size != 0)
          output << "\"" << i_vect_string->second[size - 1] << "\"";
        output << "}";
      }
    if (!output.str().empty())
      return output.str();

    if (output.str() == "")
      throw Error("LuaDefinition(string)", "Entry \"" + name
                  + "\" was not read yet in file \"" + file_path_ + "\".");

    return output.str();
  }


  //! Returns the Lua definitions of all read variables.
  /*! The variables are returned in alphabetical order.
    \return The Lua definitions of all read variables.
  */
  string Ops::LuaDefinition()
  {
    std::vector<string> name_list = GetReadEntryList();

    string output;
    std::vector<string>::iterator name;
    string previous_name = "";
    for (name = name_list.begin(); name != name_list.end(); name++)
      {
        // If '*name' is (or is part of) a new variable, inserts a new
        // line. Otherwise, '*name' is part of a table being written, so no
        // newline should be inserted.
        if (!previous_name.empty()
            && previous_name.substr(0, previous_name.find("."))
            != name->substr(0, name->find(".")))
          output += "\n";

        // Checks whether a new table is introduced. In this case, it should
        // be declared first, with "name = {}".
        if (name->find_first_of(".[") != string::npos)
          {
            string::size_type i = 0;
            string::size_type min_length
              = std::min(name->size(), previous_name.size());
            // Checks until where the previous name and the current name
            // coincide.
            while (i < min_length && previous_name[i] == (*name)[i])
              i++;
            if (i < min_length
                && (previous_name[i] == '.' || previous_name[i] == '[')
                && ((*name)[i] == '.' || (*name)[i] == '['))
              i++;
            while ((i = name->find_first_of(".[", i)) != string::npos)
              output += name->substr(0, i++) + " = {}\n";
          }

        output += LuaDefinition(*name) + "\n";

        previous_name = *name;
      }

    return output;
  }


  //! Writes the Lua definitions of all read variables.
  /*! The variables are written in alphabetical order.
    \param[in] file_name the name of the file in which the definitions are
    written.
    \note The functions are omitted.
    \warning If the output file already exists, it is cleared.
  */
  void Ops::WriteLuaDefinition(string file_name)
  {
    std::ofstream f(file_name.c_str());
    f << LuaDefinition();
    if (!f.good())
      throw Error("WriteLuaDefinition",
                  "Failed to write in \"" + file_name + "\".");
    f.close();
  }


  ///////////////////////
  // PROTECTED METHODS //
  ///////////////////////


  //! Converts an element of the stack to a reference to a single bit.
  /*! This is method is needed because a reference to an element of
    'std::vector<bool>' is not a reference to a Boolean but to a single bit.
    \param[in] index index in the stack.
    \param[out] output converted value.
    \param[in] name name of the entry.
    \return True if the conversion was successful, false otherwise.
    \note If \a name is not empty and if the conversion fails, an exception is
    raised by this method. This exception gives the name of the entry and
    states that it could not be converted. If \a name is empty, no exception
    is raised.
  */
  bool Ops::Convert(int index, std::vector<bool>::reference output,
                    string name)
  {
    if (!lua_isboolean(state_, index))
      if (name.empty())
        return false;
      else
        throw Error("Convert(vector<bool>::reference&)",
                    "The " + Entry(name) + " is not a Boolean.");

    output = static_cast<bool>(lua_toboolean(state_, index));
    return true;
  }


  //! Converts an element of the stack to a Boolean.
  /*!
    \param[in] index index in the stack.
    \param[out] output converted value.
    \param[in] name name of the entry.
    \return True if the conversion was successful, false otherwise.
    \note If \a name is not empty and if the conversion fails, an exception is
    raised by this method. This exception gives the name of the entry and
    states that it could not be converted. If \a name is empty, no exception
    is raised.
  */
  bool Ops::Convert(int index, bool& output, string name)
  {
    if (!lua_isboolean(state_, index))
      if (name.empty())
        return false;
      else
        throw Error("Convert(bool&)",
                    "The " + Entry(name) + " is not a Boolean.");

    output = static_cast<bool>(lua_toboolean(state_, index));
    return true;
  }


  //! Converts an element of the stack to an integer.
  /*!
    \param[in] index index in the stack.
    \param[out] output converted value.
    \param[in] name name of the entry.
    \return True if the conversion was successful, false otherwise.
    \note If \a name is not empty and if the conversion fails, an exception is
    raised by this method. This exception gives the name of the entry and
    states that it could not be converted. If \a name is empty, no exception
    is raised.
  */
  bool Ops::Convert(int index, int& output, string name)
  {
    if (!lua_isnumber(state_, index))
      if (name.empty())
        return false;
      else
        throw Error("Convert(int&)",
                    "The " + Entry(name) + " is not an integer.");

    double number = static_cast<double>(lua_tonumber(state_, index));
    int value = static_cast<int>(number);
    if (static_cast<double>(value) != number)
      if (name.empty())
        return false;
      else
        throw Error("Convert(int&)",
                    "The " + Entry(name) + " is not an integer.");

    output = value;
    return true;
  }


  //! Converts an element of the stack to a float.
  /*!
    \param[in] index index in the stack.
    \param[out] output converted value.
    \param[in] name name of the entry.
    \return True if the conversion was successful, false otherwise.
    \note If \a name is not empty and if the conversion fails, an exception is
    raised by this method. This exception gives the name of the entry and
    states that it could not be converted. If \a name is empty, no exception
    is raised.
  */
  bool Ops::Convert(int index, float& output, string name)
  {
    if (!lua_isnumber(state_, index))
      if (name.empty())
        return false;
      else
        throw Error("Convert(float&)",
                    "The " + Entry(name) + " is not a float.");

    output = static_cast<float>(lua_tonumber(state_, index));
    return true;
  }


  //! Converts an element of the stack to a double.
  /*!
    \param[in] index index in the stack.
    \param[out] output converted value.
    \param[in] name name of the entry.
    \return True if the conversion was successful, false otherwise.
    \note If \a name is not empty and if the conversion fails, an exception is
    raised by this method. This exception gives the name of the entry and
    states that it could not be converted. If \a name is empty, no exception
    is raised.
  */
  bool Ops::Convert(int index, double& output, string name)
  {
    if (!lua_isnumber(state_, index))
      if (name.empty())
        return false;
      else
        throw Error("Convert(double&)",
                    "The " + Entry(name) + " is not a double.");

    output = static_cast<double>(lua_tonumber(state_, index));
    return true;
  }


  //! Converts an element of the stack to a string.
  /*!
    \param[in] index index in the stack.
    \param[out] output converted value.
    \param[in] name name of the entry.
    \return True if the conversion was successful, false otherwise.
    \note If \a name is not empty and if the conversion fails, an exception is
    raised by this method. This exception gives the name of the entry and
    states that it could not be converted. If \a name is empty, no exception
    is raised.
  */
  bool Ops::Convert(int index, string& output, string name)
  {
    if (!lua_isstring(state_, index))
      if (name.empty())
        return false;
      else
        throw Error("Convert(string&)",
                    "The " + Entry(name) + " is not a string.");

    output = static_cast<string>(lua_tostring(state_, index));
    return true;
  }


  //! Retrieves a value and checks if it satisfies given constraints.
  /*! If the entry is not found, the default value is returned (if any).
    \param[in] name name of the entry.
    \param[in] constraint constraint to be satisfied.
    \param[in] default_value default value.
    \param[in] with_default is there a default value? If not, \a default_value
    is ignored.
    \param[out] value the value of the entry named \a name.
    \note The default value may not satisfy the constraint.
  */
  template<class TD, class T>
  void Ops::SetValue(string name, string constraint,
                     const TD& default_value, bool with_default,
                     T& value)
  {
    PutOnStack(Name(name));

    if (lua_isnil(state_, -1))
      if (with_default)
        {
          value = default_value;
          ClearStack();
          return;
        }
      else
        throw Error("SetValue",
                    "The " + Entry(name) + " was not found.");

    Convert(-1, value, name);

    if (!CheckConstraint(name, constraint))
      throw Error("SetValue",
                  "The " + Entry(name) + " does not satisfy "
                  + "the constraint:\n" + Constraint(constraint));

    ClearStack();

    Push(Name(name), value);
  }


  //! Retrieves a value and checks if it satisfies given constraints.
  /*! If the entry is not found, the default value is returned (if any).
    \param[in] name name of the entry.
    \param[in] constraint constraint to be satisfied.
    \param[in] default_value default value.
    \param[in] with_default is there a default value? If not, \a default_value
    is ignored.
    \param[out] value the value of the entry named \a name.
    \note The default value may not satisfy the constraint.
  */
  template<class T>
  void Ops::SetValue(string name, string constraint,
                     const std::vector<T>& default_value, bool with_default,
                     std::vector<T>& value)
  {
    PutOnStack(Name(name));

    if (lua_isnil(state_, -1))
      if (with_default)
        {
          value = default_value;
          ClearStack();
          return;
        }
      else
        throw Error("SetValue",
                    "The " + Entry(name) + " was not found.");

    if (!lua_istable(state_, -1))
      throw Error("SetValue",
                  "The " + Entry(name) + " is not a table.");

    std::vector<T> element_list;
    T element;
    std::vector<string> key_list;
    string key;
    // Now loops over all elements of the table.
    lua_pushnil(state_);
    while (lua_next(state_, -2) != 0)
      {
        // Duplicates the key and value so that 'lua_tostring' (applied to
        // them) should not interfere with 'lua_next'.
        lua_pushvalue(state_, -2);
        lua_pushvalue(state_, -2);

        if (!Convert(-2, key))
          throw Error("SetValue",
                      "Unable to read the keys of " + Entry(name) + ".");
        key_list.push_back(key);

        Convert(-1, element, name + "[" + key + "]");
        element_list.push_back(element);

        lua_pop(state_, 3);
      }

    for (int i = 0; i < int(key_list.size()); i++)
      if (!CheckConstraint(name + "[" + key_list[i] + "]", constraint))
        throw Error("SetValue",
                    "The " + Entry(name + "[" + key_list[i] + "]")
                    + " does not satisfy the constraint:\n"
                    + Constraint(constraint));

    value = element_list;

    ClearStack();

    Push(Name(name), value);
  }


  //! Prepends the prefix to an entry name.
  /*!
    \param[in] name name of the entry.
    \return The entry name with the prefix prepended.
  */
  string Ops::Name(const string& name) const
  {
    return prefix_ + name;
  }


  //! Formats the description of an entry.
  /*!
    \param[in] name name of the entry.
    \return A string with the entry name and the path to the configuration
    file, both quoted.
  */
  string Ops::Entry(const string& name) const
  {
    return "entry \"" + Name(name) + "\" in \"" + file_path_ + "\"";
  }


  //! Formats the description of a function.
  /*!
    \param[in] name name of the function.
    \return A string with the function name and the path to the configuration
    file, both quoted.
  */
  string Ops::Function(const string& name) const
  {
    return "function \"" + Name(name) + "\" in \"" + file_path_ + "\"";
  }


  //! Formats the description of a constraint.
  /*!
    \param[in] constraint the constraint to be formatted.
    \return A string with the constraint properly formatted.
  */
  string Ops::Constraint(string constraint) const
  {
    constraint = "      " + constraint;
    if (constraint.find("ops_in", 0) != string::npos)
      constraint += "\n      Note: 'ops_in(v, array)' checks whether 'v' is "
        "part of the list 'array'.";
    return constraint;
  }


  //! Iterates over the elements of \a name to put it on stack.
  /*! For instance, if "table.subtable.subsubtable[2]" is to be accessed, one
    should put "table" on top of the stack and call
    'WalkDown(subtable.subsubtable[2])'. This method will put the value of
    "table.subtable.subsubtable[2]" on top of the stack, or nil is any error
    occurred.
    \param[in] name the name of an entry that is accessible from the entry
    currently on top of the stack.
  */
  void Ops::WalkDown(string name)
  {
    if (name.empty() || lua_isnil(state_, -1))
      return;

    // The sub-entries are introduced with "." or "[i]".
    size_t end = name.find_first_of(".[");

    if (end == string::npos)
      // No more sub-entry here.
      {
        lua_pushstring(state_, name.c_str());
        lua_gettable(state_, -2);
        return;
      }

    if (name[end] == '.')
      // One step down.
      {
        lua_pushstring(state_, name.substr(0, end).c_str());
        lua_gettable(state_, -2);
        WalkDown(name.substr(end + 1).c_str());
        return;
      }

    if (name[end] == '[')
      if (end == 0)
        // Access to an element through "[i]".
        {
          // The element on stack must be a table.
          if (!lua_istable(state_, -1))
            {
              lua_pushnil(state_);
              return;
            }
          // First getting the index.
          size_t end_index = name.find_first_of("]");
          if (end_index <= end + 1 || end_index == string::npos)
            // Syntax error: "]" was not found or is misplaced.
            {
              lua_pushnil(state_);
              return;
            }
          string index_str = name.substr(end + 1, end_index - end - 1);
          // Checks whether 'index_str' is an integer.
          for (int i = 0; i < int(index_str.size()); i++)
            if (!isdigit(index_str[i]))
              {
                lua_pushnil(state_);
                return;
              }
          std::istringstream str(index_str);
          int index;
          str >> index;
          // Now getting the element of index 'index'.
          lua_rawgeti(state_, -1, index);
          // And preparing for the next step down.
          string next_name = name.substr(end_index + 1).c_str();
          if (!next_name.empty() && next_name[0] == '.')
            // Removes the dot in the first place.
            next_name = next_name.substr(1);
          WalkDown(next_name);
          return;
        }
      else
        // One step down before addressing "[i]" in the next call to
        // 'WalkDown'.
        {
          lua_pushstring(state_, name.substr(0, end).c_str());
          lua_gettable(state_, -2);
          if (name[end] == '.')
            WalkDown(name.substr(end + 1).c_str());
          else
            WalkDown(name.substr(end).c_str());
          return;
        }
  }


  //! Checks whether \a name is of type 'T'.
  /*! On exit, the value of the entry (if it exists) is on the stack.
    \param[in] name the name of the entry whose type is checked.
    \param[in] value anything: it is used to determine the type.
    \return True if the entry is of type 'T', false otherwise.
    \note The prefix is prepended to \a name. If \a name does not exist, an
    exception is raised.
  */
  template<class T>
  bool Ops::IsParam(string name, T& value)
  {
    PutOnStack(Name(name));

    if (lua_isnil(state_, -1))
      throw Error("Is(string)",
                  "The " + Entry(name) + " was not found.");

    return Convert(-1, value);
  }


  //! Checks whether \a name is a table of 'T'.
  /*!
    \param[in] name the name of the entry whose type is checked.
    \param[in] value anything: it is used to determine the type.
    \return True if the entry is a table of 'T', false otherwise.
    \note The prefix is prepended to \a name. If \a name does not exist, an
    exception is raised.
  */
  template<class T>
  bool Ops::IsParam(string name, std::vector<T>& value)
  {
    PutOnStack(Name(name));

    if (lua_isnil(state_, -1))
      throw Error("IsParam",
                  "The " + Entry(name) + " was not found.");

    if (!lua_istable(state_, -1))
      return false;

    T element;
    // Now loops over all elements of the table.
    lua_pushnil(state_);
    while (lua_next(state_, -2) != 0)
      {
        // Duplicates the value so that 'lua_tostring' (applied to them)
        // should not interfere with 'lua_next'.
        lua_pushvalue(state_, -1);

        if (!Convert(-1, element))
          return false;
        lua_pop(state_, 2);
      }

    return true;
  }


  //! Stores the value of an entry.
  /*!
    \param[in] name the name of the entry.
    \param[in] value the value of the entry.
  */
  void Ops::Push(string name, const bool& value)
  {
    read_bool[name] = value;
  }


  //! Stores the value of an entry.
  /*!
    \param[in] name the name of the entry.
    \param[in] value the value of the entry.
  */
  void Ops::Push(string name, const int& value)
  {
    read_int[name] = value;
  }


  //! Stores the value of an entry.
  /*!
    \param[in] name the name of the entry.
    \param[in] value the value of the entry.
  */
  void Ops::Push(string name, const float& value)
  {
    read_float[name] = value;
  }


  //! Stores the value of an entry.
  /*!
    \param[in] name the name of the entry.
    \param[in] value the value of the entry.
  */
  void Ops::Push(string name, const double& value)
  {
    read_double[name] = value;
  }


  //! Stores the value of an entry.
  /*!
    \param[in] name the name of the entry.
    \param[in] value the value of the entry.
  */
  void Ops::Push(string name, const string& value)
  {
    read_string[name] = value;
  }


  //! Stores the value of an entry.
  /*!
    \param[in] name the name of the entry.
    \param[in] value the value of the entry.
  */
  void Ops::Push(string name, const std::vector<bool>& value)
  {
    read_vect_bool[name] = value;
  }


  //! Stores the value of an entry.
  /*!
    \param[in] name the name of the entry.
    \param[in] value the value of the entry.
  */
  void Ops::Push(string name, const std::vector<int>& value)
  {
    read_vect_int[name] = value;
  }


  //! Stores the value of an entry.
  /*!
    \param[in] name the name of the entry.
    \param[in] value the value of the entry.
  */
  void Ops::Push(string name, const std::vector<float>& value)
  {
    read_vect_float[name] = value;
  }


  //! Stores the value of an entry.
  /*!
    \param[in] name the name of the entry.
    \param[in] value the value of the entry.
  */
  void Ops::Push(string name, const std::vector<double>& value)
  {
    read_vect_double[name] = value;
  }


  //! Stores the value of an entry.
  /*!
    \param[in] name the name of the entry.
    \param[in] value the value of the entry.
  */
  void Ops::Push(string name, const std::vector<string>& value)
  {
    read_vect_string[name] = value;
  }


  //! Pushes all keys of a map into a vector.
  /*!
    \param[in] input the map whose keys should be pushed.
    \param[in,out] vect the vector to which the map keys are pushed back.
  */
  template<class TK, class T>
  void Ops::AppendKey(const std::map<TK, T>& input, std::vector<TK>& vect)
  {
    typename std::map<TK, T>::const_iterator i;
    for (i = input.begin(); i != input.end(); i++)
      vect.push_back(i->first);
  }


}


#define OPS_FILE_CLASSOPS_CXX
#endif
