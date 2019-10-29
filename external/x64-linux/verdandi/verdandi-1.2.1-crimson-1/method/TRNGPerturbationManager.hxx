// Copyright (C) 2010-2011 INRIA
// Author(s): Kévin Charpentier, Vivien Mallet, Anne Tilloy
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


#ifndef VERDANDI_FILE_METHOD_TRNGPERTURBATIONMANAGER_HXX
#define VERDANDI_FILE_METHOD_TRNGPERTURBATIONMANAGER_HXX


#include "BasePerturbationManager.hxx"
#include "trng/yarn5.hpp"


namespace Verdandi
{


    /////////////////////////////
    // TRNGPERTURBATIONMANAGER //
    /////////////////////////////


    //! This class generates random samples using TRNG.
    class TRNGPerturbationManager:
        public BasePerturbationManager<TRNGPerturbationManager>
    {
    protected:
        //! Normal random number generator.
        trng::yarn5* nrng_;

    public:

        /*** Constructors and destructor ***/

        TRNGPerturbationManager();
        TRNGPerturbationManager(string configuration_file);
        ~TRNGPerturbationManager();

        /*** Methods ***/

        void Initialize(string configuration_file);
        void Initialize(VerdandiOps& configuration_stream);
        void Finalize();

        double Normal(double mean, double variance,
                      Vector<double, VectFull>& parameter);
        double LogNormal(double mean, double variance,
                         Vector<double, VectFull>& parameter);

        template <class T0, class T1,
                  class Prop0, class Allocator0>
        void Normal(Matrix<T0, Prop0, RowSymPacked, Allocator0> variance,
                    Vector<double, VectFull>& parameter,
                    Vector<T1, VectFull, Allocator0>& sample);

        template <class T0, class Prop0, class Allocator0,
                  class T1, class Allocator1>
        void LogNormal(Matrix<T0, Prop0, RowSymPacked, Allocator0> variance,
                       Vector<double, VectFull>& parameter,
                       Vector<T1, VectFull, Allocator1>& output);

        template <class T0,
                  class T1, class Allocator1>
        void NormalHomogeneous(T0 variance,
                               Vector<double, VectFull>& parameter,
                               Vector<T1, VectFull, Allocator1>& output);

        template <class T0,
                  class T1, class Allocator1>
        void LogNormalHomogeneous(T0 variance,
                                  Vector<double, VectFull>& parameter,
                                  Vector<T1, VectFull, Allocator1>& output);
        template <class T0,
                  class T1, class Allocator0>
        bool NormalClipping(Vector<T0, VectFull>& diagonal,
                            Vector<double, VectFull>& parameter,
                            Vector<T1, VectFull, Allocator0>& output);
    };


} // namespace Verdandi.


#endif
