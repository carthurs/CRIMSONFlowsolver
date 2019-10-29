// Copyright (C) 2008-2010, INRIA
// Author(s): Vivien Mallet, Claire Mouton, Marc Fragu
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


#ifndef VERDANDI_FILE_PYTHON_VERDANDI_CPP


#include "verdandi.def"

#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK
#define SELDON_WITH_SUPERLU

#define VERDANDI_DENSE

#include "Verdandi.hxx"
#include "VerdandiBase.cxx"

#include "seldon/SeldonSolver.hxx"

#include "model/QuadraticModel.cxx"
#include "model/ClampedBar.cxx"
#include "model/PythonModel.cxx"
#include "observation_manager/GridToNetworkObservationManager.cxx"
#include "observation_manager/LinearObservationManager.cxx"
#include "observation_manager/PythonObservationManager.cxx"
#include "method/OptimalInterpolation.cxx"
#include "method/ForwardDriver.cxx"
#include "method/ReducedOrderExtendedKalmanFilter.cxx"
#include "share/OutputSaver.cxx"

namespace Verdandi
{


    template class VSWIG_MODEL;
    template class VSWIG_MODEL1;
    class VSWIG_MODEL2;
    template class VSWIG_GRID_TO_NETWORK_OBSERVATION;
    template class VSWIG_LINEAR_OBSERVATION;
    class VSWIG_LINEAR_OBSERVATION2;
    template class VSWIG_METHOD;
    template class VSWIG_METHOD1;
    template class VSWIG_METHOD2;
    template class VSWIG_METHOD3;
    template class VSWIG_METHOD4;
    template class VSWIG_METHOD5;


}

#define VERDANDI_FILE_PYTHON_VERDANDI_CPP
#endif
