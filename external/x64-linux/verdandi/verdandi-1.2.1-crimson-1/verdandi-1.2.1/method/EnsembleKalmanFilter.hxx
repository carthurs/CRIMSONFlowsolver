// Copyright (C) 2011 INRIA
// Author(s): Kévin Charpentier, Vivien Mallet
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

#ifndef VERDANDI_FILE_METHOD_ENSEMBLEKALMANFILTER_HXX


namespace Verdandi
{


    //////////////////////////
    // ENSEMBLEKALMANFILTER //
    //////////////////////////


    //! This class implements the ensemble Kalman filter.
    template <class T, class Model, class ObservationManager,
              class PerturbationManager>
    class EnsembleKalmanFilter: public VerdandiBase
    {
    public:
        //! Type of the model state vector.
        typedef typename Model::state model_state;
        //! Type of the observation vector.
        typedef typename ObservationManager::observation
        observation;
        //! Type of an uncertain parameter.
        typedef typename Model::uncertain_parameter uncertain_parameter;
        //! Type of the ensemble of state vectors.
        typedef vector<model_state> ensemble;

    protected:

        /*** Main components ***/

        //! Underlying model.
        Model model_;

        //! Observation manager.
        ObservationManager observation_manager_;

        //! Pertubation managers.
        PerturbationManager perturbation_manager_;

        //! The number of members in the ensemble.
        int Nmember_;

        //! Local number of ensemble members (for parallel computing).
        int Nlocal_member_;

        //! The ensemble state vectors.
        ensemble ensemble_;

        //! Ensemble for the parameters.
        vector<vector<uncertain_parameter> > parameter_;

        //! Number of parameters to be perturbed.
        int Nparameter_;

        /*** Configuration ***/

        //! Path to the configuration file.
        string configuration_file_;
        //! Path to the model configuration file.
        string model_configuration_file_;
        //! Path to the configuration file for the perturbation manager.
        string perturbation_manager_configuration_file_;
        //! Path to the configuration file for the observation manager.
        string observation_configuration_file_;

        //! Display options.
        map<string, bool> option_display_;

        //! Dimension of the reduced state.
        int Nstate_;
        //! Dimension of the full state.
        int Nfull_state_;
        //! Number of observations.
        int Nobservation_;
        //! Should an analysis be computed at the first step?
        bool analyze_first_step_;
        //! Current time.
        double time_;

#if defined(VERDANDI_WITH_MPI)

        /*** Parallel settings ***/

        //! Process rank.
        int rank_;
        //! Number of processes.
        int Nprocess_;
#endif

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

    public:

        /*** Constructor and destructor ***/

        EnsembleKalmanFilter();
        ~EnsembleKalmanFilter();

        /*** Methods ***/

        void Initialize(string configuration_file,
                        bool initialize_model = true,
                        bool initialize_observation_manager = true,
                        bool initialize_perturbation_manager = true);
        void Initialize(VerdandiOps& configuration,
                        bool initialize_model = true,
                        bool initialize_observation_manager = true,
                        bool initialize_perturbation_manager = true);

        void InitializeEnsemble();

        void InitializeStep();

        void Forward();
        void Analyze();

        void FinalizeStep();
        void Finalize();

        bool HasFinished() const;

        // Access methods.
        Model& GetModel();
        ObservationManager& GetObservationManager();
        PerturbationManager& GetPerturbationManager();
        OutputSaver& GetOutputSaver();
        string GetName() const;

        void Message(string message);

    protected:

        template <class T0, class Allocator0>
        void Fill(Vector<T0, Collection, Allocator0>& in, string pdf);
        template <class T0, class Storage0, class Allocator0>
        void Fill(Vector<T0, Storage0, Allocator0>& in, string pdf);
        template <class T0, class Storage0, class Allocator0>
        void SetDimension(Vector<T0, Storage0, Allocator0>& in,
                          Vector<T0, Storage0, Allocator0>& out);
        template <class T0, class Allocator0>
        void SetDimension(Vector<T0, Collection, Allocator0>& in,
                          Vector<T0, Collection, Allocator0>& out);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_ENSEMBLEKALMANFILTER_HXX
#endif
