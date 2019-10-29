// Copyright (C) 2009 INRIA
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


#ifndef VERDANDI_FILE_METHOD_FORWARDDRIVER_HXX


namespace Verdandi
{


    ///////////////////
    // FORWARDDRIVER //
    ///////////////////


    //! This class simply performs a forward simulation.
    template <class Model>
    class ForwardDriver: public VerdandiBase
    {

    public:
        //! Type of the model state vector.
        typedef typename Model::state model_state;

    protected:

        //! Underlying model.
        Model model_;

        //! Iteration.
        int iteration_;
        //! Time vector.
        Vector<double> time_;

        /*** Configuration ***/

        //! Path to the configuration file.
        string configuration_file_;
        //! Path to the model configuration file.
        string model_configuration_file_;

        //! Should the iterations be displayed?
        bool show_iteration_;
        //! Should the current time be displayed?
        bool show_time_;

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

    public:

        /*** Constructor and destructor ***/

        ForwardDriver();
        ~ForwardDriver();

        /*** Methods ***/

        void Initialize(string configuration_file,
                        bool initialize_model = true);
        void Initialize(VerdandiOps& configuration,
                        bool initialize_model = true);
        void InitializeStep();
        void Forward();
        void FinalizeStep();
        void Finalize();

        bool HasFinished();

        // Access methods.
        Model& GetModel();
        OutputSaver& GetOutputSaver();
        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_FORWARDDRIVER_HXX
#endif
