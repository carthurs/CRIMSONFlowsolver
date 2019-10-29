// Copyright (C) 2008-2009 INRIA
// Author(s): Vivien Mallet, Claire Mouton
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


#ifndef VERDANDI_FILE_VERDANDIHEADER_HXX


//! The namespace of the data assimilation library Verdandi.
namespace Verdandi
{


} // namespace Verdandi.


#include <iostream>
#include <fstream>
#include <map>
#include <list>
#include <cmath>
#include <utility>


//////////////////
// DEBUG LEVELS //
//////////////////


#ifdef VERDANDI_DEBUG_LEVEL_4
#ifndef VERDANDI_DEBUG_LEVEL_3
#define VERDANDI_DEBUG_LEVEL_3
#endif
#define SELDON_DEBUG_LEVEL_4
#endif

#ifdef VERDANDI_DEBUG_LEVEL_3
#ifndef VERDANDI_DEBUG_LEVEL_2
#define VERDANDI_DEBUG_LEVEL_2
#endif
#define SELDON_DEBUG_LEVEL_3
#endif

#ifdef VERDANDI_DEBUG_LEVEL_2
#ifndef VERDANDI_DEBUG_LEVEL_1
#define VERDANDI_DEBUG_LEVEL_1
#endif
#ifndef VERDANDI_CHECK_DIMENSIONS
#define VERDANDI_CHECK_DIMENSIONS
#endif
#define SELDON_DEBUG_LEVEL_2
#endif

#ifdef VERDANDI_DEBUG_LEVEL_1
#ifndef VERDANDI_DEBUG_LEVEL_0
#define VERDANDI_DEBUG_LEVEL_0
#endif
#define SELDON_DEBUG_LEVEL_1
#ifndef VERDANDI_CHECK_IO
#define VERDANDI_CHECK_IO
#endif
#endif

#ifdef VERDANDI_DEBUG_LEVEL_0
#define SELDON_DEBUG_LEVEL_0
#ifndef VERDANDI_DEBUG_LEVEL_2
#define VERDANDI_IS_ACTIVE false
#endif
#endif


#ifdef VERDANDI_WITH_ABORT
#define SELDON_WITH_ABORT
#endif

#ifndef VERDANDI_WITH_ABORT
#define OPS_WITH_EXCEPTION
#endif


// Convenient macros to catch exceptions.
#ifndef VERDANDI_TRY
#define VERDANDI_TRY try {
#endif
#ifndef VERDANDI_END
#define VERDANDI_END                                            \
  }                                                             \
    catch(Verdandi::Error& Err)                                 \
      {                                                         \
        Err.CoutWhat();                                         \
        return 1;                                               \
      }                                                         \
    catch(Seldon::Error& Err)                                   \
      {                                                         \
        Err.CoutWhat();                                         \
        return 1;                                               \
      }                                                         \
    catch(::Ops::Error& Err)                                    \
      {                                                         \
        Err.CoutWhat();                                         \
        return 1;                                               \
      }                                                         \
    catch(std::exception& Err)                                  \
      {                                                         \
        std::cout << "C++ exception: "                          \
            << Err.what() << std::endl;                         \
        return 1;                                               \
      }                                                         \
    catch(std::string& str)                                     \
      {                                                         \
        std::cout << str << std::endl;                          \
        return 1;                                               \
      }                                                         \
    catch(const char* str)                                      \
      {                                                         \
        std::cout << str << std::endl;                          \
        return 1;                                               \
      }                                                         \
    catch(...)                                                  \
      {                                                         \
        std::cout << "Unknown exception..." << std::endl;       \
        return 1;                                               \
      }
#endif


#include "seldon/SeldonHeader.hxx"
#include "seldon/vector/Vector2.hxx"
#include "seldon/vector/Vector3.hxx"


namespace Verdandi
{


    using namespace std;

    using namespace Seldon;

    using Seldon::to_num;
    using Seldon::to_str;


} // namespace Verdandi.


#ifndef VERDANDI_PYTHON_VERSION
#define VERDANDI_PYTHON_VERSION 2.6
#endif

#include "share/VerdandiOps.hxx"
#include "share/Logger.hxx"
#include "share/Error.hxx"
#include "share/UsefulFunction.hxx"
#include "share/MessageHandler.hxx"
#include "share/VerdandiBase.hxx"
#include "share/OutputSaver.hxx"

#ifdef VERDANDI_SPARSE
#define VERDANDI_TANGENT_LINEAR_OPERATOR_SPARSE
#define VERDANDI_OBSERVATION_ERROR_SPARSE
#define VERDANDI_STATE_ERROR_SPARSE
#endif

#ifdef VERDANDI_DENSE
#define VERDANDI_TANGENT_LINEAR_OPERATOR_DENSE
#define VERDANDI_STATE_ERROR_DENSE
#endif

#define VERDANDI_FILE_VERDANDIHEADER_HXX
#endif
