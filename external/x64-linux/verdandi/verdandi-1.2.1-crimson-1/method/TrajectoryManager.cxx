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
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_FILE_METHOD_TRAJECTORYMANAGER_CXX

#include "TrajectoryManager.hxx"

#include "seldon/vector/VectorCollection.cxx"


namespace Verdandi
{


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////


    //! Main constructor.
    template <class T, class Model>
    TrajectoryManager<T, Model>::TrajectoryManager()
    {
    }


    //! Destructor.
    template <class T, class Model>
    TrajectoryManager<T, Model>::~TrajectoryManager()
    {
        Deallocate();
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the trajectory manager.
    /*
      \param[in] checkpoint_recording_mode checkpoint recording mode
      (memory or disk).
      \param[in] checkpoint_recording_file checkpoint recording file.
      \param[in] loaded_trajectory_recording_mode loaded trajectory
      recording mode.
      \param[in] loaded_trajectory_recording_file loaded trajectory
      recording file.
      \param[in] Nskip_checkpoint recording period for checkpoint.
    */
    template <class T, class Model>
    void TrajectoryManager<T, Model>
    ::Initialize(string checkpoint_recording_mode,
                 string checkpoint_recording_file,
                 string loaded_trajectory_recording_mode,
                 string loaded_trajectory_recording_file,
                 int Nskip_checkpoint)
    {
        checkpoint_recording_mode_ = checkpoint_recording_mode;
        checkpoint_recording_file_ = checkpoint_recording_file;
        loaded_trajectory_recording_mode_ = loaded_trajectory_recording_mode;
        loaded_trajectory_recording_file_ = loaded_trajectory_recording_file;
        Nskip_checkpoint_ = Nskip_checkpoint;
        Nsave_call_ = 0;
        Ncheckpoint_ = 0;
        Nstate_ = 0;
        checkpoint_index_ = -1;
        loaded_trajectory_index_ = -1;
        if (checkpoint_recording_mode_ == "disk")
            EmptyFile(checkpoint_recording_file);
        if (loaded_trajectory_recording_mode == "disk")
            EmptyFile(loaded_trajectory_recording_mode);
    }


    //! It appends \a state at a given time \a time to the checkpoints.
    /*!
      \param state a given state.
      \param time the given time.
    */
    template <class T, class Model>
    void TrajectoryManager<T, Model>::Save(model_state& state,
                                           double time)
    {
        if ((Nsave_call_ % Nskip_checkpoint_) == 0)
        {
            if (checkpoint_recording_mode_ == "memory")
                checkpoint_.push_back(state);
            else
            {
                fstream file_stream;
                file_stream.open(checkpoint_recording_file_.c_str());
#ifdef VERDANDI_CHECK_IO
                // Checks if the file was opened.
                if (!file_stream.is_open())
                    throw ErrorIO("TrajectoryManager<T, Model>"
                                  "::Save", "Unable to open file \""
                                  + checkpoint_recording_file_ + "\".");
#endif
                streampos
                    position = Ncheckpoint_ *
                    (state.GetM() * sizeof(T) + sizeof(int));
                file_stream.seekg(position);
                state.Write(file_stream, true);
            }
            checkpoint_time_.PushBack(time);
            Ncheckpoint_++;
        }
        Nsave_call_++;
    }


    //! Sets the time of states to be loaded.
    /*!
      \param[in] model the model.
      \param[in] time a given time.
    */
    template <class T, class Model>
    void TrajectoryManager<T, Model>::SetTime(Model& model, double time)
    {
        Nstate_ = model.GetNstate();
        if (checkpoint_index_ >= 0
            && loaded_trajectory_index_ >= 0
            && loaded_time_(loaded_trajectory_index_) == time)
            return;

        if (checkpoint_index_ >= 0
            && loaded_trajectory_index_ > 0
            && loaded_time_(loaded_trajectory_index_ - 1) == time)
        {
            loaded_trajectory_index_--;
            return;
        }

        if (checkpoint_index_ >= 0
            && checkpoint_time_(checkpoint_index_) <= time
            && (checkpoint_index_ == Ncheckpoint_ - 1 ||
                time < checkpoint_time_(checkpoint_index_ + 1)))
        {
            bool available_time;
            for (int i = 0; i < loaded_time_.GetM(); i++)
            {
                available_time = false;
                if (loaded_time_(i) == time)
                {
                    available_time = true;
                    loaded_trajectory_index_ = i;
                    break;
                }
            }
            if (!available_time)
                throw ErrorProcessing("void TrajectoryManager<T, Model>"
                                      "::SetTime(Model& model, double time)",
                                      "No state was saved at time: "
                                      + to_str(time) + ".");
            return;
        }
        if (checkpoint_index_ > 0
            && checkpoint_time_(checkpoint_index_ - 1) <= time
            && time < checkpoint_time_(checkpoint_index_))
            checkpoint_index_--;
        else
        {
            bool available_time = false;
            for (int i = 0; i < Ncheckpoint_; i++)
            {
                if (checkpoint_time_(i) <= time &&
                    (i == Ncheckpoint_ - 1 ||
                     time < checkpoint_time_(i + 1)))

                {
                    available_time = true;
                    checkpoint_index_ = i;
                    break;
                }
            }
            if (!available_time)
                throw ErrorProcessing("void TrajectoryManager<T, Model>"
                                      "::SetTime(Model& model, double time)",
                                      "No state was saved at time: "
                                      + to_str(time) + ".");
        }

        model_state saved_state(model.GetNstate());
        double saved_time;
        model_state& tmp = model.GetState();
        Copy(tmp, saved_state);
        saved_time = model.GetTime();

        if (loaded_trajectory_recording_mode_ == "memory")
            loaded_trajectory_.clear();
        else
            EmptyFile(loaded_trajectory_recording_file_);
        loaded_time_.Clear();
        model.SetTime(checkpoint_time_(checkpoint_index_));
        if (checkpoint_recording_mode_ == "memory")
        {
            Copy(checkpoint_[checkpoint_index_], tmp);
            model.StateUpdated();
        }
        else
        {
            ifstream file_stream;
            file_stream.open(checkpoint_recording_file_.c_str());
#ifdef VERDANDI_CHECK_IO
            // Checks if the file was opened.
            if (!file_stream.is_open())
                throw ErrorIO("TrajectoryManager<T, Model>"
                              "::SetTime", "Unable to open file \""
                              + checkpoint_recording_file_ + "\".");
#endif
            model_state input_data;
            streampos position;
            position = checkpoint_index_ *
                (model.GetNstate() * sizeof(T) + sizeof(int));
            file_stream.seekg(position);
            input_data.Read(file_stream);
            file_stream.close();
            Copy(input_data, tmp);
            model.StateUpdated();
        }

        if (loaded_trajectory_recording_mode_ == "memory")
        {
            double tmp_time;
            for (int t = 0; t < Nskip_checkpoint_ && !model.HasFinished();
                 t++)
            {
                model_state& state = model.GetState();
                tmp_time = model.GetTime();
                loaded_trajectory_.push_back(state);
                loaded_time_.PushBack(tmp_time);
                model.Forward();
                if (tmp_time == time)
                    loaded_trajectory_index_ = t;
            }
        }
        else
        {
            fstream file_stream;
            file_stream.open(loaded_trajectory_recording_file_.c_str());
#ifdef VERDANDI_CHECK_IO
            // Checks if the file was opened.
            if (!file_stream.is_open())
                throw ErrorIO("TrajectoryManager<T, Model>"
                              "::SetTime", "Unable to open file \""
                              + loaded_trajectory_recording_file_ + "\".");
#endif
            double tmp_time;
            for (int t = 0; t < Nskip_checkpoint_ && !model.HasFinished();
                 t++)
            {
                model_state& state = model.GetState();
                tmp_time = model.GetTime();
                state.Write(file_stream, true);
                loaded_time_.PushBack(tmp_time);
                model.Forward();
                if (tmp_time == time)
                    loaded_trajectory_index_ = t;
            }
            file_stream.close();
        }
        Copy(saved_state, tmp);
        model.StateUpdated();
        model.SetTime(saved_time);
    }


    //! Returns the state at current time.
    /*!
      \return The current state.
    */
    template <class T, class Model>
    typename TrajectoryManager<T, Model>::model_state&
    TrajectoryManager<T, Model>::GetState()
    {
        if (loaded_trajectory_recording_mode_ == "memory")
            return loaded_trajectory_[loaded_trajectory_index_];

        if (Nstate_ == 0)
            throw ErrorProcessing("TrajectoryManager<T, Model>::GetState()",
                                  "Nstate = 0");
        ifstream file_stream;
        file_stream.open(loaded_trajectory_recording_file_.c_str());
#ifdef VERDANDI_CHECK_IO
        // Checks if the file was opened.
        if (!file_stream.is_open())
            throw ErrorIO("TrajectoryManager<T, Model>"
                          "::GetState", "Unable to open file \""
                          + loaded_trajectory_recording_file_ + "\".");
#endif
        streampos position;
        position = loaded_trajectory_index_ *
            (Nstate_ * sizeof(T) + sizeof(int));
        file_stream.seekg(position);
        input_data_.Read(file_stream);
        file_stream.close();
        return input_data_;
    }


    //! Returns the current time.
    /*!
      \return The current time.
    */
    template <class T, class Model>
    double TrajectoryManager<T, Model>::GetTime()
    {
        return loaded_time_(loaded_trajectory_index_);
    }


    //! Clears inner vector collection.
    template <class T, class Model>
    void TrajectoryManager<T, Model>::Deallocate()
    {
        checkpoint_.clear();
        loaded_trajectory_.clear();
        checkpoint_time_.Clear();
        loaded_time_.Clear();
        Ncheckpoint_ = 0;
        checkpoint_index_ = -1;
        loaded_trajectory_index_ = -1;
        Nsave_call_ = 0;
        if (checkpoint_recording_mode_ == "disk")
            EmptyFile(checkpoint_recording_file_);
        if (loaded_trajectory_recording_mode_ == "disk")
            EmptyFile(loaded_trajectory_recording_file_);
    }


    //! Empties the recording file.
    /*!
      \param[in] file_name path to the file to be emptied.
    */
    template <class T, class Model>
    void TrajectoryManager<T, Model>::EmptyFile(string file_name)
    {
        ofstream file(file_name.c_str());
        if (file)
            file.close();
        else
            throw ErrorIO("TrajectoryManager::EmptyFile(string)",
                          "Cannot open file \"" + file_name + "\"." );
    }



} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_TRAJECTORYMANAGER_CXX
#endif
