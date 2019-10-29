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


%module ops
%{
#include "OpsHeader.hxx"
  %}

%include "typemaps.i"
%include "std_string.i"
%include "std_vector.i"
using namespace std;

namespace std
{
  %template(VectBool) vector<bool>;
  %template(VectInt) vector<int>;
  %template(VectFloat) vector<float>;
  %template(VectDouble) vector<double>;
  %template(VectString) vector<string>;
}

%include "Error.hxx"
%exception
{
  try
    {
      $action
	}
  catch(Ops::Error& e)
    {
      PyErr_SetString(PyExc_Exception, e.What().c_str());
      return NULL;
    }
  catch(std::exception& e)
    {
      PyErr_SetString(PyExc_Exception, e.what());
      return NULL;
    }
  catch(std::string& s)
    {
      PyErr_SetString(PyExc_Exception, s.c_str());
      return NULL;
    }
  catch(const char* s)
    {
      PyErr_SetString(PyExc_Exception, s);
      return NULL;
    }
  catch(...)
    {
      PyErr_SetString(PyExc_Exception, "Unknown exception...");
      return NULL;
    }
}

%include "OpsHeader.hxx"
%include "ClassOps.hxx"


%define OPS_INSTANTIATE_ELEMENT(suffix, type)
%template(Get ## suffix) Get<type >;
%template(Apply ## suffix) Apply<type >;
%template(Is ## suffix) Is<type >;
%enddef

%define OPS_INSTANTIATE_CROSSED_ELEMENT(suffix0, type0, suffix1, type1)
%template(Apply ## suffix0 ## suffix1) Apply<type0, type1 >;
%enddef

%define OPS_INSTANTIATE_VECTOR(suffix, type)
%template(Get ## suffix) Get<type >;
%template(Is ## suffix) Is<type >;
%enddef

namespace Ops
{

  %extend Ops
  {
    OPS_INSTANTIATE_ELEMENT(Bool, bool);
    OPS_INSTANTIATE_ELEMENT(Int, int);
    OPS_INSTANTIATE_ELEMENT(Float, float);
    OPS_INSTANTIATE_ELEMENT(Double, double);
    OPS_INSTANTIATE_ELEMENT(String, string);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Bool, bool, Bool, bool);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Int, int, Bool, bool);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Float, float, Bool, bool);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Double, double, Bool, bool);
    OPS_INSTANTIATE_CROSSED_ELEMENT(String, string, Bool, bool);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Bool, bool, Int, int);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Int, int, Int, int);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Float, float, Int, int);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Double, double, Int, int);
    OPS_INSTANTIATE_CROSSED_ELEMENT(String, string, Int, int);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Bool, bool, Float, float);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Int, int, Float, float);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Float, float, Float, float);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Double, double, Float, float);
    OPS_INSTANTIATE_CROSSED_ELEMENT(String, string, Float, float);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Bool, bool, Double, double);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Int, int, Double, double);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Float, float, Double, double);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Double, double, Double, double);
    OPS_INSTANTIATE_CROSSED_ELEMENT(String, string, Double, double);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Bool, bool, String, string);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Int, int, String, string);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Float, float, String, string);
    OPS_INSTANTIATE_CROSSED_ELEMENT(Double, double, String, string);
    OPS_INSTANTIATE_CROSSED_ELEMENT(String, string, String, string);
    OPS_INSTANTIATE_VECTOR(VectBool, std::vector<bool>);
    OPS_INSTANTIATE_VECTOR(VectInt, std::vector<int>);
    OPS_INSTANTIATE_VECTOR(VectFloat, std::vector<float>);
    OPS_INSTANTIATE_VECTOR(VectDouble, std::vector<double>);
    OPS_INSTANTIATE_VECTOR(VectString, std::vector<string>);
  };

}
