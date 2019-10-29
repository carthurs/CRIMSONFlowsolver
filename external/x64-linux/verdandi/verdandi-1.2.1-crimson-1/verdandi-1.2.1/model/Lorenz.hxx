// Copyright (C) 2008-2010 INRIA
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


#ifndef VERDANDI_FILE_MODEL_LORENZ_HXX


namespace Verdandi
{


    //////////////////
    // LORENZ MODEL //
    //////////////////


    //! This class is a Lorenz model.
    /*!
      \tparam T the type of floating-point numbers.
    */
    template <class T>
    class Lorenz: public VerdandiBase
    {
    public:
        typedef T value_type;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T& reference;
        typedef const T& const_reference;
        typedef Matrix<T> state_error_variance;
        typedef Vector<T> state_error_variance_row;
        typedef Vector<T> state;
        typedef Matrix<T> matrix_state_observation;

    protected:

        //! Prandtl number.
        T Pr_;
        //! Rayleigh number.
        T Ra_;
        //! Third parameter.
        T b_;

        //! First variable of the system.
        T X_;
        //! Second variable of the system.
        T Y_;
        //! Third variable of the system.
        T Z_;

        //! Backup of the first variable of the system.
        T X_tmp_;
        //! Backup of the second variable of the system.
        T Y_tmp_;

        //! Time step.
        double Delta_t_;

        //! Final time of the simulation.
        double final_time_;

        //! Current time.
        double time_;

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

        //! State vector.
        state state_;

    public:
        // Constructor and destructor.
        Lorenz();
        Lorenz(string configuration_file);
        ~Lorenz();
        // Initializations.
        void Initialize(string configuration_file);
        void InitializeStep();

        // Processing.
        void Forward();
        bool HasFinished() const;
        void Save();

        void FinalizeStep();
        void Finalize();

        // Access methods.
        T GetX() const;
        T GetY() const;
        T GetZ() const;
        T GetDelta_t() const;
        double GetTime() const;
        void SetTime(double time);
        int GetNstate() const;
        state& GetState();
        void StateUpdated();
        state& GetFullState();
        void FullStateUpdated();
        void GetStateErrorVarianceRow(int row, state_error_variance_row&
                                      error_covariance_row);
        const state_error_variance& GetStateErrorVariance() const;
        bool IsErrorSparse() const;

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_LORENZ_HXX
#endif
