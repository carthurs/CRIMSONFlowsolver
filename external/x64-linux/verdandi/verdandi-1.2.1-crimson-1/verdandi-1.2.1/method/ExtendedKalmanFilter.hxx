// Copyright (C) 2008-2009 INRIA
// Author(s): Marc Fragu, Vivien Mallet, Claire Mouton
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


#ifndef VERDANDI_FILE_EXTENDEDKALMANFILTER_HXX


namespace Verdandi
{


    //////////////////////////
    // EXTENDEDKALMANFILTER //
    //////////////////////////


    //! This class implements the extended Kalman filter.
    template <class T, class Model, class ObservationManager>
    class ExtendedKalmanFilter: public VerdandiBase
    {

    public:
        //! Type of a row of the background error variance.
        typedef typename Model::state_error_variance_row
        model_state_error_variance_row;
        //! Type of the model state vector.
        typedef typename Model::state model_state;
        //! Type of the model/observation crossed matrix.
        typedef typename Model::matrix_state_observation
        matrix_state_observation;
        //! Type of the background error variance.
        typedef typename Model::state_error_variance
        model_state_error_variance;
        //! Type of the tangent linear model.
        typedef typename Model::tangent_linear_operator
        model_tangent_linear_operator;
        //! Type of the tangent linear observation operator.
        typedef typename ObservationManager
        ::tangent_linear_operator observation_tangent_linear_operator;
        //! Type of a row of the tangent linear observation operator.
        typedef typename ObservationManager::tangent_linear_operator_row
        observation_tangent_linear_operator_row;
        //! Type of the observation vector.
        typedef typename ObservationManager::observation
        observation;

    protected:

        /*** Main components ***/

        //! Underlying model.
        Model model_;

        //! Observation manager.
        ObservationManager observation_manager_;

        //! Background error covariance matrix (B).
        model_state_error_variance state_error_variance_;

        /*** Configuration ***/

        //! Path to the configuration file.
        string configuration_file_;
        //! Path to the model configuration file.
        string model_configuration_file_;
        //! Path to the configuration file for the observation manager.
        string observation_configuration_file_;

        //! Display options.
        map<string, bool> option_display_;

        //! Dimension of the state.
        int Nstate_;
        //! Number of observations.
        int Nobservation_;
        //! Should an analysis be computed at the first step?
        bool analyze_first_step_;
        //! Current time.
        double time_;

        //! Computation mode for BLUE: "vector" or "matrix".
        string blue_computation_;
        //! Computation mode for covariance: "vector" or "matrix".
        string covariance_computation_;

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

    public:

        /*** Constructors and destructor ***/

        ExtendedKalmanFilter();
        ~ExtendedKalmanFilter();

        /*** Methods ***/

        void Initialize(string configuration_file,
                        bool initialize_model = true,
                        bool initialize_observation_manager = true);
        void Initialize(VerdandiOps& configuration,
                        bool initialize_model = true,
                        bool initialize_observation_manager = true);

        void InitializeStep();

        void Forward();
        void Analyze();

        void FinalizeStep();
        void Finalize();

        void PropagateCovarianceMatrix();
        void PropagateCovarianceMatrix_vector();
        void PropagateCovarianceMatrix_matrix();

        void ComputeBLUE(const observation& innovation, model_state& state);

        bool HasFinished();

        // Access methods.
        Model& GetModel();
        ObservationManager& GetObservationManager();
        OutputSaver& GetOutputSaver();

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_EXTENDEDKALMANFILTER_HXX
#endif
