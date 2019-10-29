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


// This file is used to instantiate explicitely the template methods, so as to
// build a compiled library and the Python interface.


#ifndef OPS_FILE_OPS_CPP

#include "Ops.hxx"

#define OPS_INSTANTIATE_ELEMENT(type)                           \
  template type Ops::Get(string);                               \
  template type Ops::Get(string, string);                       \
  template type Ops::Get(string, string, const type&);          \
  template type Ops::Apply(string name, const type& arg0);      \
  template type Ops::Apply(string name, const type& arg0,       \
                           const type& arg1);                   \
  template type Ops::Apply(string name, const type& arg0,       \
                           const type& arg1, const type& arg2); \
  template type Ops::Apply(string name, const type& arg0,       \
                           const type& arg1, const type& arg2,  \
                           const type& arg3);                   \
  template type Ops::Apply(string name, const type& arg0,       \
                           const type& arg1, const type& arg2,  \
                           const type& arg3, const type& arg4); \
  template bool Ops::Is<type >(string);                         \

#define OPS_INSTANTIATE_CROSSED_ELEMENT(type0, type1)           \
  template void Ops::Apply(string name,                         \
                           const std::vector<type0>& in,        \
                           std::vector<type1>& out);            \

#define OPS_INSTANTIATE_VECTOR(type)                    \
  template type Ops::Get(string);                       \
  template type Ops::Get(string, string);               \
  template type Ops::Get(string, string, const type&);  \
  template bool Ops::Is<type >(string);                 \

namespace Ops
{
  OPS_INSTANTIATE_ELEMENT(bool);
  OPS_INSTANTIATE_ELEMENT(int);
  OPS_INSTANTIATE_ELEMENT(float);
  OPS_INSTANTIATE_ELEMENT(double);
  OPS_INSTANTIATE_ELEMENT(string);
  OPS_INSTANTIATE_CROSSED_ELEMENT(bool, bool);
  OPS_INSTANTIATE_CROSSED_ELEMENT(int, bool);
  OPS_INSTANTIATE_CROSSED_ELEMENT(float, bool);
  OPS_INSTANTIATE_CROSSED_ELEMENT(double, bool);
  OPS_INSTANTIATE_CROSSED_ELEMENT(string, bool);
  OPS_INSTANTIATE_CROSSED_ELEMENT(bool, int);
  OPS_INSTANTIATE_CROSSED_ELEMENT(int, int);
  OPS_INSTANTIATE_CROSSED_ELEMENT(float, int);
  OPS_INSTANTIATE_CROSSED_ELEMENT(double, int);
  OPS_INSTANTIATE_CROSSED_ELEMENT(string, int);
  OPS_INSTANTIATE_CROSSED_ELEMENT(bool, float);
  OPS_INSTANTIATE_CROSSED_ELEMENT(int, float);
  OPS_INSTANTIATE_CROSSED_ELEMENT(float, float);
  OPS_INSTANTIATE_CROSSED_ELEMENT(double, float);
  OPS_INSTANTIATE_CROSSED_ELEMENT(string, float);
  OPS_INSTANTIATE_CROSSED_ELEMENT(bool, double);
  OPS_INSTANTIATE_CROSSED_ELEMENT(int, double);
  OPS_INSTANTIATE_CROSSED_ELEMENT(float, double);
  OPS_INSTANTIATE_CROSSED_ELEMENT(double, double);
  OPS_INSTANTIATE_CROSSED_ELEMENT(string, double);
  OPS_INSTANTIATE_CROSSED_ELEMENT(bool, string);
  OPS_INSTANTIATE_CROSSED_ELEMENT(int, string);
  OPS_INSTANTIATE_CROSSED_ELEMENT(float, string);
  OPS_INSTANTIATE_CROSSED_ELEMENT(double, string);
  OPS_INSTANTIATE_CROSSED_ELEMENT(string, string);
  OPS_INSTANTIATE_VECTOR(std::vector<bool>);
  OPS_INSTANTIATE_VECTOR(std::vector<int>);
  OPS_INSTANTIATE_VECTOR(std::vector<float>);
  OPS_INSTANTIATE_VECTOR(std::vector<double>);
  OPS_INSTANTIATE_VECTOR(std::vector<string>);
}


#define OPS_FILE_OPS_CPP
#endif
