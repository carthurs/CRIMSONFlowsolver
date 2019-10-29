// Copyright (C) 2009 INRIA
// Author(s): Vivien Mallet, Claire Mouton, Anne Tilloy
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

#ifndef VERDANDI_FILE_METHOD_MONTECARLO_HXX


namespace Verdandi
{


    ////////////////
    // MONTECARLO //
    ////////////////


    //! This class performs allows to perform Monte Carlo simulations.
    /*! The class performs a single simulation with perturbed data. In order
      to complete a full Monte Carlo simulation, one has to launch several
      simulations using this class.
    */
    template <class T, class ClassModel, class PerturbationManager>
    class MonteCarlo: public VerdandiBase
    {

    public:
        typedef typename ClassModel::state model_state;
        typedef typename ClassModel::uncertain_parameter uncertain_parameter;

    protected:

        //! Underlying model.
        ClassModel model_;

        //! Pertubation managers.
        PerturbationManager perturbation_manager_;

        //! Perturbations vectors.
        vector<uncertain_parameter> perturbation_;

        //! Iteration.
        int iteration_;
        //! Time vector.
        Vector<double> time_;

        /*** Configuration ***/

        //! Path to the configuration file.
        string configuration_file_;
        //! Path to the model configuration file.
        string model_configuration_file_;
        //! Path to the configuration file for the perturbation manager.
        string perturbation_manager_configuration_file_;

        //! Should the iterations be displayed?
        bool show_iteration_;
        //! Should the current time be displayed?
        bool show_time_;

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

    public:

        /*** Constructor and destructor ***/

        MonteCarlo();
        ~MonteCarlo();
        template <class T0, class Storage0, class Allocator0>
        void Clear(Vector<T0, Storage0, Allocator0>& V);
        template <class T0, class Allocator0>
        void Clear(Vector<T0, Collection, Allocator0>& V);

        /*** Main methods ***/

        template <class T0, class Storage0, class Allocator0>
        void SetDimension(Vector<T0, Storage0, Allocator0>& in,
                          Vector<T0, Storage0, Allocator0>& out);
        template <class T0, class Allocator0>
        void SetDimension(Vector<T0, Collection, Allocator0>& in,
                          Vector<T0, Collection, Allocator0>& out);

        template <class T0, class Allocator0>
        void Fill(Vector<T0, Collection, Allocator0>& in, string pdf);
        template <class T0, class Storage0, class Allocator0>
        void Fill(Vector<T0, Storage0, Allocator0>& in, string pdf);

        void Initialize(string configuration_file,
                        bool initialize_model = true,
                        bool initialize_perturbation_manager = true);
        void Initialize(VerdandiOps& configuration,
                        bool initialize_model = true,
                        bool initialize_perturbation_manager = true);
        void InitializeStep();
        void Forward();
        void FinalizeStep();
        void Finalize();

        bool HasFinished();

        /*** Access methods ***/

        ClassModel& GetModel();
        OutputSaver& GetOutputSaver();
        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_MONTECARLO_HXX
#endif
