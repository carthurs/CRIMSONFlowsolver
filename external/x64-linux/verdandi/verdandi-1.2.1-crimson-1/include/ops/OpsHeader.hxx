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


#ifndef OPS_FILE_OPSHEADER_HXX


//! Ops namespace.
namespace Ops
{
}


extern "C"
{
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

#ifndef OPS_WITH_EXCEPTION
#define OPS_WITH_ABORT
#endif

#ifndef DISP
#define DISP(x) std::cout << #x ": " << x << std::endl
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <map>

using std::string;

#include "Error.hxx"
#include "ClassOps.hxx"


#define OPS_FILE_OPSHEADER_HXX
#endif
