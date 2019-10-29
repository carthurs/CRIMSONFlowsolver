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


#ifndef OPS_FILE_CLASSOPS_HXX


namespace Ops
{

  class Ops
  {
  protected:
    //! Path to the configuration file.
    string file_path_;
    //! Lua state.
    lua_State* state_;
    //! Prefix to be prepended to the entries names.
    string prefix_;

    //! Names and values of all Booleans read in the file.
    std::map<string, bool> read_bool;
    //! Names and values of all integer read in the file.
    std::map<string, int> read_int;
    //! Names and values of all floats read in the file.
    std::map<string, float> read_float;
    //! Names and values of all doubles read in the file.
    std::map<string, double> read_double;
    //! Names and values of all strings read in the file.
    std::map<string, string> read_string;
    //! Names and values of all vectors of Booleans read in the file.
    std::map<string, std::vector<bool> > read_vect_bool;
    //! Names and values of all vectors of integers read in the file.
    std::map<string, std::vector<int> > read_vect_int;
    //! Names and values of all vectors of floats read in the file.
    std::map<string, std::vector<float> > read_vect_float;
    //! Names and values of all vectors of doubles read in the file.
    std::map<string, std::vector<double> > read_vect_double;
    //! Names and values of all vectors of strings read in the file.
    std::map<string, std::vector<string> > read_vect_string;

  public:
    // Constructor and destructor.
    Ops();
    Ops(string file_path);
    ~Ops();

    // Main methods.
    void Open(string file_path, bool close_state = true);
    void Reload(bool close_state = true);
    void Close();
    template<class TD, class T>
    void
    Set(string name, string constraint, const TD& default_value, T& value);
    template<class T>
    void Set(string name, string constraint, T& value);
    template<class T>
    void Set(string name, T& value);
    template<class T>
    T Get(string name);
    template<class T>
    T Get(string name, string constraint);
    template<class T>
    T Get(string name, string constraint, const T& default_value);
    template<class Tin, class Tout>
    void Apply(string name, const std::vector<Tin>& in,
               std::vector<Tout>& OUTPUT);
    template<class T>
    T Apply(string name, const T& arg0);
    template<class T>
    T Apply(string name, const T& arg0, const T& arg1);
    template<class T>
    T Apply(string name, const T& arg0, const T& arg1, const T& arg2);
    template<class T>
    T Apply(string name, const T& arg0, const T& arg1, const T& arg2,
            const T& arg3);
    template<class T>
    T Apply(string name, const T& arg0, const T& arg1, const T& arg2,
            const T& arg3, const T& arg4);
    std::vector<string> GetEntryList(string name = "");
    bool CheckConstraint(string name, string constraint);
    bool CheckConstraintOnValue(string value, string constraint);
    void PutOnStack(string name);
    bool Exists(string name);
    void PushOnStack(bool value);
    void PushOnStack(int value);
    void PushOnStack(float value);
    void PushOnStack(double value);
    void PushOnStack(string value);
    template<class T>
    void PushOnStack(const std::vector<T>& v);
    template<class T>
    bool Is(string name);
    bool IsTable(string name);
    bool IsFunction(string name);
    void ClearStack();

    void DoFile(string file_path);
    void DoString(string expression);

    // Access methods.
    string GetFilePath() const;
    lua_State* GetState();
#ifndef SWIG
    const lua_State* GetState() const;
#endif
    string GetPrefix() const;
    void SetPrefix(string prefix);
    void ClearPrefix();
    std::vector<string> GetReadEntryList();
    void UpdateLuaDefinition();
    string LuaDefinition(string name);
    string LuaDefinition();
    void WriteLuaDefinition(string file_name);

  protected:
    bool Convert(int index, std::vector<bool>::reference output,
                 string name = "");
    bool Convert(int index, bool& output, string name = "");
    bool Convert(int index, int& output, string name = "");
    bool Convert(int index, float& output, string name = "");
    bool Convert(int index, double& output, string name = "");
    bool Convert(int index, string& output, string name = "");
    template<class TD, class T>
    void SetValue(string name, string constraint,
                  const TD& default_value, bool with_default, T& value);
    template<class T>
    void SetValue(string name, string constraint,
                  const std::vector<T>& default_value, bool with_default,
                  std::vector<T>& value);
    string Constraint(string constraint) const;
    string Name(const string& name) const;
    string Entry(const string& name) const;
    string Function(const string& name) const;
    void WalkDown(string name);
    template<class T>
    bool IsParam(string name, T& value);
    template<class T>
    bool IsParam(string name, std::vector<T>& value);
    void Push(string name, const bool& value);
    void Push(string name, const int& value);
    void Push(string name, const float& value);
    void Push(string name, const double& value);
    void Push(string name, const string& value);
    void Push(string name, const std::vector<bool>& value);
    void Push(string name, const std::vector<int>& value);
    void Push(string name, const std::vector<float>& value);
    void Push(string name, const std::vector<double>& value);
    void Push(string name, const std::vector<string>& value);
    template<class TK, class T>
    void AppendKey(const std::map<TK, T>& input, std::vector<TK>& vect);
  };

}


#define OPS_FILE_CLASSOPS_HXX
#endif
