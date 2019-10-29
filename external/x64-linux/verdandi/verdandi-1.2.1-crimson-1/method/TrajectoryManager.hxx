// Copyright (C) 2011 INRIA
// Author(s): Marc Fragu
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


#ifndef VERDANDI_FILE_METHOD_TRAJECTORYMANAGER_HXX


namespace Verdandi
{


    ///////////////////////
    // TRAJECTORYMANAGER //
    ///////////////////////


    //! This class manages model trajectory.
    template <class T, class Model>
    class TrajectoryManager
    {

    public:
        //! Type of the model state vector.
        typedef typename Model::state model_state;

    protected:

        //! Checkpoint recording mode.
        string checkpoint_recording_mode_;
        //! Checkpoint recording file name.
        string checkpoint_recording_file_;
        //! Loaded trajectory recording mode.
        string loaded_trajectory_recording_mode_;
        //! Loaded trajectory aving file name.
        string loaded_trajectory_recording_file_;
        //! Number of call to save method.
        int Nsave_call_;
        //! Model Nstate.
        int Nstate_;
        //! Recording period.
        int Nskip_checkpoint_;
        //! Number of recording in first forward loop.
        int Ncheckpoint_;
        //! Trajectory saved in first forward loop.
        vector<model_state> checkpoint_;
        //! Times associated with \a checkpoint_ .
        Vector<double> checkpoint_time_;

        //! Trajectory saved in backward loop.
        vector<model_state> loaded_trajectory_;
        //! Times associated with \a loaded_trajectory_ .
        Vector<double> loaded_time_;
        //! Checkpoint index.
        int checkpoint_index_;
        //! State loaded from file.
        model_state input_data_;
        //! Loaded trajectory index.
        int loaded_trajectory_index_;

    public:

        /*** Constructor and destructor ***/

        TrajectoryManager();
        ~TrajectoryManager();

        /*** Methods ***/

        void Initialize(string checkpoint_recording_mode = "memory",
                        string checkpoint_recording_file = "",
                        string loaded_trajectory_recording_mode = "memory",
                        string loaded_trajectory_recording_file = "",
                        int Nskip_checkpoint = 1);
        void Save(model_state& state, double time);
        void SetTime(Model& model, double time);

        model_state& GetState();
        double GetTime();
        void Deallocate();
        void EmptyFile(string file_name);

    };


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_TRAJECTORYMANAGER_HXX
#endif
