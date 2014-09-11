#ifndef VERDANDI_FILE_MODEL_SimvascularVerdandiModel_CXX
#define VERDANDI_FILE_MODEL_SimvascularVerdandiModel_CXX

#include "SimvascularVerdandiModel.hxx"

//! Checks for existance of file
/*!
      \param[in] const std::string& name

      This function simply takes a file name
      and attempts to open it.
      If successful, it returns true,
      otherwise, it returns false.
 */
inline bool exists_test1 (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

namespace Verdandi {


////////////////////////////////
// CONSTRUCTOR AND DESTRUCTOR //
////////////////////////////////


/*
 *    This is the default constructor for SimvascularVerdandiModel.
 *    Currently, it sets initial values for the
 *    data members.
 */
SimvascularVerdandiModel::SimvascularVerdandiModel()
:   gat(NULL),
 	Nreduced_(0),
 	Nstate_(0),
 	Nstate_local_(0),
 	state_reduced_start_local_(-1),
 	shared_parts_size_(0),
 	nreduced_has_wall_parameters_(0),
 	nreduced_has_coupled_parameters_(0),
 	cp_rcr_estimate_resistance_(0),
 	cp_rcr_estimate_compliance_(0),
 	cp_rcr_estimate_prox_resistance_(0),
 	cp_rcr_estimate_pout_(0),
 	cp_rcr_estimate_pstates_(0),
 	rank_(0),
 	numProcs_(1),
 	iNewComm_C_(MPI_COMM_WORLD)
{
}


/*
 *    This is the default destructor for SimvascularVerdandiModel.
 */
SimvascularVerdandiModel::~SimvascularVerdandiModel() {
	Eoutfile_.close();
	//MPI_Finalize();
}


////////////////
// INITIALIZE //
////////////////



/*!
      \param[in] configuration_file configuration file.

      This function opens the Verdandi configuration file
      and reads in configuration settings for the model.
      It then calls the default Initialize.
 */
void SimvascularVerdandiModel::Initialize(string configuration_file) {

	VerdandiOps configuration(configuration_file);

	configuration.SetPrefix("simvascular_model.");
	configuration.Set("state_reduced_has_wall_parameters",nreduced_has_wall_parameters_);
	configuration.Set("state_reduced_has_coupled_parameters",nreduced_has_coupled_parameters_);

	/* 
	   Note that in the *.lua file, the RCR_parameters_info is a block
       i.e. RCR_parameters_info = {estimate_resistance, ..., ...}
	*/

	if (nreduced_has_coupled_parameters_) {
		configuration.Set("RCR_parameters_info.estimate_resistance",cp_rcr_estimate_resistance_);
		configuration.Set("RCR_parameters_info.estimate_compliance",cp_rcr_estimate_compliance_);
		configuration.Set("RCR_parameters_info.estimate_prox_resistance",cp_rcr_estimate_prox_resistance_);
		configuration.Set("RCR_parameters_info.estimate_pout",cp_rcr_estimate_pout_);
		configuration.Set("RCR_parameters_info.estimate_pstates",cp_rcr_estimate_pstates_);

		configuration.Set("RCR_parameters_info.resistance_included",cp_rcr_include_resistance_);
		configuration.Set("RCR_parameters_info.compliance_included",cp_rcr_include_compliance_);
		configuration.Set("RCR_parameters_info.prox_resistance_included",cp_rcr_include_prox_resistance_);

        configuration.Set("Heart_parameters_info.estimate_emax",cp_hrt_estimate_emax_);		
        configuration.Set("Heart_parameters_info.estimate_tmax",cp_hrt_estimate_tmax_);		
        configuration.Set("Heart_parameters_info.estimate_trel",cp_hrt_estimate_trel_);		                
        configuration.Set("Heart_parameters_info.estimate_vlv",cp_hrt_estimate_vlv_);		
	}

	configuration.Set("error_statistics.state_error_variance",state_error_variance_value_);

    Initialize();

}


/*!
      This is the function that sets up
      the MPI communicator,
      calls the PHASTA routines for file input,
      partitioning, and initialization.
      This is also where the vectors
      of augmented state components are populated
      and the sizes of the augmented state,
      local part of the augmented state,
      and the estimated part of the augmented state
      are tabulated.
 */
void SimvascularVerdandiModel::Initialize() {

	char pathToProcsCaseDir[100];

	// Get pointer to the single instance of SimvascularGlobalArrayTransfer
	gat = SimvascularGlobalArrayTransfer::Get();

	// save the communicator
	iNewComm_C_ = MPI_COMM_WORLD;
	newcom.iNewComm = MPI_Comm_c2f(iNewComm_C_); // modifies newcom in fortran common block

	MPI_Comm_size(iNewComm_C_, &numProcs_);
	MPI_Comm_rank(iNewComm_C_, &rank_);

	// read configuration file
	input_fform();

	// set iestimator to 1, i.e. estimator 
	// int exists in the common block
	nomodule.iestimator = int(1);

	// Preprocess data and run the problem
	// Partition the problem to the correct number of processors
	if( rank_ == 0 )
	{
		//std::cout << "number of procs per simulation " << numprocs_perparticle << std::endl;
		Partition_Problem( numProcs_ );
	}

	MPI_Barrier(iNewComm_C_);

	sprintf(pathToProcsCaseDir,"%d-procs-case",numProcs_);
	chdir(pathToProcsCaseDir);

	input(&numProcs_, &rank_);
    proces();
    itrdrv_init(); // initialize solver

    // organize the arrays from PHASTA into the SimvascularAugStatePart class
    BuildAugmentedState();

    int estimated_total = 0, state_total = 0;

    std::cout << "at rank " << rank_ << " number of states " << dstrb_parts_.size() << std::endl;

    for(auto it = dstrb_parts_.begin(); it != dstrb_parts_.end(); ++it) {
    	std::cout << "at rank " << rank_ << " " << it->getName().c_str() << " " << it->getSize() <<  std::endl;
    	state_total += (int)it->getSize();
    }

    if (rank_ == numProcs_ - 1) {

    	for (auto it = shared_parts_.begin(); it != shared_parts_.end(); ++it) {
    		std::cout << "at rank " << rank_ << " " << it->getName().c_str() << " " << it->getSize() << std::endl;

    		int first_estimated = it->getFirstEstimated();

    		if (first_estimated >= 0 && state_reduced_start_local_ < 0)
    			state_reduced_start_local_ = (int)state_total+first_estimated;

    		estimated_total += (int)it->getNumEstimated();
    		state_total += (int)it->getSize();
    	}
    }

    for (auto it = shared_parts_.begin(); it != shared_parts_.end(); ++it)
    	shared_parts_size_ += it->getSize();

    Nstate_local_ = state_total; // the size of the on-proc duplicated_state_

    Nreduced_ = estimated_total; // the size of the estimated portion

    // we want Nreduced_ on every proc
    MPI_Allreduce(&Nreduced_, &Nreduced_, 1, MPI_INT, MPI_SUM, iNewComm_C_);

    std::cout << "at rank " << rank_ << " # of distributed states: " << Nstate_local_ << std::endl;
    std::cout << "at rank " << rank_ << " # of shared states: " << shared_parts_size_ << std::endl;
    std::cout << "at rank " << rank_ << " # of estimated states: " << Nreduced_ << std::endl;
    std::cout << "at rank " << rank_ << " state_reduced_start " << state_reduced_start_local_ << std::endl;

	// Compute the global state size
	MPI_Allreduce(&Nstate_local_, &Nstate_, 1, MPI_INT, MPI_SUM, iNewComm_C_);




	// Make sure the state array is distributed with the correct number
	// of elements per processor that is given by Nstate_local_
	duplicated_state_.Reallocate(Nstate_, Nstate_local_);


	// ***
	// *** set the indices in duplicated_state_ to get/set
	// ***

	int state_start, state_end, icounter;

	duplicated_state_.GetProcessorRange(state_start, state_end);

	icounter = state_start;

	for(auto it = dstrb_parts_.begin(); it != dstrb_parts_.end(); ++it)
		for (std::size_t jj = 0; jj < it->getSize(); jj++)
			it->addDuplicatedStateIndex(icounter++);

	if (rank_ == numProcs_ - 1)
		for (auto it = shared_parts_.begin(); it != shared_parts_.end(); ++it)
			for (std::size_t jj = 0; jj < it->getSize(); jj++)
				it->addDuplicatedStateIndex(icounter++);


    // file output for initial step if starting from scratch
    std::string Efilename = "estimated_reduced_state.dat";

    if (exists_test1(Efilename)) { // still need to finish this
        // The file exists, and is open for input
    	std::cout << "we can start from the previous estimates" << std::endl;
    } else {
    	std::cout << "we need to write the new file" << std::endl;
    }

    if (rank_ == numProcs_ - 1) {

    	Eoutfile_.open(Efilename.c_str(),std::ofstream::app);
    	Eoutfile_ << std::scientific << std::setprecision( std::numeric_limits<double>::digits10 );

    }

}


/*!
 *    The purpose of this function is to
 *    build up the augmented state.
 *    Specifically, it populates
 *    the dstrb_parts_ and shared_parts_ vectors
 *    with SimvascularAugStatePart instances,
 *    corresponding to various PHASTA Fortran arrays.
 */
void SimvascularVerdandiModel::BuildAugmentedState() {

	SimvascularAugStatePart state_part;


	// velocity and pressure field
	// here we are loop through the nodes only on this processor, i.e. unique/not shared
	// the number of these nodes is nshguniq, this was set-up in partition KDL NAN
	state_part.Initialize("solution");
	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {
		int actualIdx = (gat->pointerMapInt_["local index of unique nodes"])[unitIdx];
		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond
			state_part.addDataPointer(&gat->pointerMapDP_["solution"][varIdx * conpar.nshg + actualIdx-1]);
			state_part.addIsEstimated(0);
		}
	}
	dstrb_parts_.push_back(state_part);
	state_part.Clear();

	// acceleration fields
	state_part.Initialize("time derivative of solution");
	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {
		int actualIdx = (gat->pointerMapInt_["local index of unique nodes"])[unitIdx];
		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond
			state_part.addDataPointer(&gat->pointerMapDP_["time derivative of solution"][varIdx * conpar.nshg + actualIdx-1]);
			state_part.addIsEstimated(0);
		}
	}
	dstrb_parts_.push_back(state_part);
	state_part.Clear();

	if (nomodule.ideformwall > 0) {

		// displacement field
		state_part.Initialize("displacement");
		for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {
			int actualIdx = (gat->pointerMapInt_["local index of unique nodes"])[unitIdx];
			for(int varIdx=0; varIdx < 3; varIdx++) {
				state_part.addDataPointer(&gat->pointerMapDP_["displacement"][varIdx * conpar.nshg + actualIdx-1]);
				state_part.addIsEstimated(0);
			}
		}
		dstrb_parts_.push_back(state_part);
		state_part.Clear();

	}

	// this is where the estimated parts of the augmented state will go

	// lumped parameter state (for RCR, just P_current)
	if (grcrbccom.numGRCRSrfs > 0) {
		state_part.Initialize("RCR P_current");
		for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++) {
			state_part.addDataPointer(&gat->pointerMapDP_["WindkesselRCR_P"][surfIdx]);
			if (cp_rcr_estimate_pstates_)
				state_part.addIsEstimated(1);
			else
				state_part.addIsEstimated(0);
		}
		shared_parts_.push_back(state_part);
		state_part.Clear();
	}


	// lumped parameter state (for heart model, just VLV)
	if (nomodule.iheart > 0) {
		state_part.Initialize("Heart LV Volume");		
		state_part.addDataPointer(&gat->pointerMapDP_["Heart_LV_Vol"][0]);
        std::cout << state_part.getName() << " " << gat->pointerMapDP_["Heart_LV_Vol"][0] << endl;

		// estimating this variable
		if (cp_hrt_estimate_vlv_){
			state_part.addIsEstimated(1);
		} else {
			state_part.addIsEstimated(0);
		}
		shared_parts_.push_back(state_part);
		state_part.Clear();
	}

    // ***
    // *** now come the parameters, above are the variables required to restart the simulation, i.e. P^{t_{n}}, etc.
    // *** below are the estimated parameters, these have to be collected together
    // ***

	// for vessel wall regional properties (E)
	if (nreduced_has_wall_parameters_ && nomodule.ideformwall > 0 && nomodule.numWallRegions > 0) {
		state_part.Initialize("E");
		for(int parIdx = 0; parIdx < nomodule.numWallRegions; parIdx++) {
			state_part.addDataPointer(&gat->pointerMapDP_["Vessel Wall Young's Modulus"][parIdx]);
			state_part.addIsEstimated(1);
		}
		shared_parts_.push_back(state_part);
		state_part.Clear();
	}

	// for lumped parameter model parameters
	// for RCR
	if (nreduced_has_coupled_parameters_) {

		// for RCR C
		if (cp_rcr_estimate_compliance_) {

			state_part.Initialize("RCR C");
			for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++) {

				if (cp_rcr_include_compliance_[parIdx]) {
					state_part.addDataPointer(&gat->pointerMapDP_["WindkesselRCR_Params"][parIdx*3+1]);
					//std::cout << gat->pointerMapDP_["WindkesselRCR_Params"][parIdx*3+1] << endl;
					state_part.addIsEstimated(1);
				}

			}
			shared_parts_.push_back(state_part);
			state_part.Clear();

		}

		// for RCR R2
		if (cp_rcr_estimate_resistance_) {

			state_part.Initialize("RCR R2");
			for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++) {

				if (cp_rcr_include_resistance_[parIdx]) {
					state_part.addDataPointer(&gat->pointerMapDP_["WindkesselRCR_Params"][parIdx*3+2]);
					//std::cout << gat->pointerMapDP_["WindkesselRCR_Params"][parIdx*3+2] << endl;
					state_part.addIsEstimated(1);
				}

			}
			shared_parts_.push_back(state_part);
			state_part.Clear();

		}

		// for RCR R1
		if (cp_rcr_estimate_prox_resistance_) {

			state_part.Initialize("RCR R1");
			for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++) {

				if (cp_rcr_include_prox_resistance_[parIdx]) {
					state_part.addDataPointer(&gat->pointerMapDP_["WindkesselRCR_Params"][parIdx*3+0]);
					//std::cout << gat->pointerMapDP_["WindkesselRCR_Params"][parIdx*3+0] << endl;
					state_part.addIsEstimated(1);
				}

			}
			shared_parts_.push_back(state_part);
			state_part.Clear();

		}

		// for RCR Pdist
		if (cp_rcr_estimate_pout_) {

			state_part.Initialize("RCR Pdist");

			state_part.addDataPointer(&gat->pointerMapDP_["WindkesselRCR_Pdist"][0]);
			std::cout << gat->pointerMapDP_["WindkesselRCR_Pdist"][0] << endl;
			state_part.addIsEstimated(1);

			shared_parts_.push_back(state_part);
			state_part.Clear();

		}

        // for heart model EMax 
		if (cp_hrt_estimate_emax_) {

			// give it a name
			state_part.Initialize("Heart EMax");

            // single parameter, will be in the 0 index and print out 
			state_part.addDataPointer(&gat->pointerMapDP_["Heart_EMax"][0]);           
			std::cout << state_part.getName() << " " << gat->pointerMapDP_["Heart_EMax"][0] << endl;
			
			// set this parameter to be estimated and add to shared_parts vector
			state_part.addIsEstimated(1);
			shared_parts_.push_back(state_part);

			// clear object
			state_part.Clear();
		}

        // for heart model TMax 
		if (cp_hrt_estimate_tmax_) {

			// give it a name
			state_part.Initialize("Heart TMax");

            // single parameter, will be in the 0 index and print out 
			state_part.addDataPointer(&gat->pointerMapDP_["Heart_TMax"][0]);           
			std::cout << state_part.getName() << " " << gat->pointerMapDP_["Heart_TMax"][0] << endl;
			
			// set this parameter to be estimated and add to shared_parts vector
			state_part.addIsEstimated(1);
			shared_parts_.push_back(state_part);

			// clear object
			state_part.Clear();
		}

        // for heart model TRel 
		if (cp_hrt_estimate_trel_) {

			// give it a name
			state_part.Initialize("Heart TRel");

            // single parameter, will be in the 0 index and print out 
			state_part.addDataPointer(&gat->pointerMapDP_["Heart_TRel"][0]);           
			std::cout << state_part.getName() << " " << gat->pointerMapDP_["Heart_TRel"][0] << endl;
			
			// set this parameter to be estimated and add to shared_parts vector
			state_part.addIsEstimated(1);
			shared_parts_.push_back(state_part);

			// clear object
			state_part.Clear();
		}

	}


	for (auto it = dstrb_parts_.begin(); it != dstrb_parts_.end(); ++it)
		dstrb_parts_map_[it->getName()] = &(*it);

	for (auto it = shared_parts_.begin(); it != shared_parts_.end(); ++it)
		shared_parts_map_[it->getName()] = &(*it);

}

/*
 *    Things that need to be done
 *    to initalize the first time step
 *    of the model should go here.
 */
void SimvascularVerdandiModel::InitializeFirstStep() {

}


/*
 *    Things that need to be done
 *    to initialize a time step
 *    shoudl go here.
 */
void SimvascularVerdandiModel::InitializeStep() {

}


/*
 *    Things that need to be done
 *    to finalizes the model.
 *    Currently, this calls the finalize
 *    portion of the itrdrv Fortran routine
 */
void SimvascularVerdandiModel::Finalize() {

	itrdrv_finalize();

}


////////////////
// PROCESSING //
////////////////



/*!
 *    Propagate the model forward by one time step
 *
 *    \f[X_{k+1} = A_k(X_k, \theta_k)\,.\f]
 *
 */
void SimvascularVerdandiModel::Forward() {

	itrdrv_iter_init();

	itrdrv_iter_step();

	itrdrv_iter_finalize();

}


/*
 *    This simply calls the PHASTA Fortran
 *    routine itrdrv_iter_finalize
 *    to allow moving on to the
 *    next time step
 */
void SimvascularVerdandiModel::FinalizeStep() {

	itrdrv_iter_finalize();

}




/*!
      \return True if the simulation is done, false otherwise.

      This simply checks if the time step
      is greater or equal than the
      number of time steps specified in solver.inp
      Previously the name loop control
      was done entirely in itrdrv.
 */
bool SimvascularVerdandiModel::HasFinished() const {

	return timdat.istep >= inpdat.nstep[0];
}


///////////////
// OPERATORS //
///////////////



/*! The current state of the model is modified.
      \param[in] x a vector.
      \param[in] forward Boolean to indicate if the model has to go on to the
      next step.
      \param[in] preserve_state Boolean to indicate if the model state has to
      be preserved.

      This function takes an arbitrary input state,
      copies the relevant variables to
      the appropriate PHASTA Fortran arrays,
      calls the PHASTA functions to apply
      the forward dynamics, and then copies
      the resulting updated state
      back into the input state.
      Unlike the Forward function,
      ApplyOperator does not call
      itrdrv_iter_finalize to finalize
      the time step.  This is because multiple calls
      to ApplyOperator will happen (i.e., multiple particles)
      before we will finalize the time step.
 */
void SimvascularVerdandiModel::ApplyOperator(state& x,
		bool forward, bool preserve_state) {

	//double saved_time = 0;
	//state saved_state;
	//if (!forward)
	//  saved_time = GetTime();
    //if (preserve_state)
    //  GetStateCopy(saved_state);

	duplicated_state_.Copy(x); // copies x into the duplicated state vector
	StateUpdated();            // updates actual model state with duplicated state vector

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
	{
		stringstream s_temp;
		s_temp << this->rank_;
		string outname = "x_internal_before_apply-"+s_temp.str()+".dat";
		ofstream outfile(outname.c_str(), ofstream::app);
		outfile.precision(std::numeric_limits<double>::digits10);

		string outnameb = "x_internal_before_apply-"+s_temp.str()+".bin";
		ofstream outfileb(outnameb.c_str(), ofstream::app);
		outfileb.precision(15);

		Vector <double> local_X;
		double* local_temp;
		int nxtemp = duplicated_state_.GetLocalM();
		VecGetArray(duplicated_state_.GetPetscVector(), &local_temp);  // get data pointer from petsc vector
		local_X.SetData(nxtemp,local_temp);  // set pointer for local vector

		local_X.WriteText(outfile);
		local_X.Write(outfileb);
		outfile << std::endl;
		local_X.Nullify();
	}
#endif

	itrdrv_iter_init();  // advances the model forward

	itrdrv_iter_step();  // note that ForwardFinalize is not called here

	x.Copy(GetState());        // copies the actual model state (via duplicated state) into x

}


////////////////////
// ACCESS METHODS //
////////////////////



/*!
      \return: The current time in the model.
 */
double SimvascularVerdandiModel::GetTime() const {

	return (double)(timdat.lstep);
}



/*!
      \param[in] time a given time.
 */
void SimvascularVerdandiModel::SetTime(double time) {
	throw ErrorUndefined("SimvascularVerdandiModel"
			"::SetTime(double time)");
}



/*!
      \return The augmented state vector size.
 */
int SimvascularVerdandiModel::GetNstate() const {
	return Nstate_;
}



/*!
      \return The size of the local augmented state vector.
 */
int SimvascularVerdandiModel::GetLocalNstate() const {
	return Nstate_local_;
}


/*!
      \return The starting index of the local reduced state vector.
 */
int SimvascularVerdandiModel::GetLocalReducedStart() const {
	return state_reduced_start_local_;
}



/*!
      \return: a reference to "duplicated state" vector

      This function first updates
      the internal duplicated_state_ with the values
      from the appropriate PHASTA Fortran arrays,
      via the pointers stored in
      the dstrb_parts_ and shared_parts_
      vectors.  Note that this function
      is completely agnostic to the entries in these
      vectors, which were populated
      in BuildAugmentedState.
      It then returns a reference to the
      update duplicated_state_.
 */
SimvascularVerdandiModel::state& SimvascularVerdandiModel::GetState() {

	if (rank_ == 0)
		std::cout << "getting state ";

	for(auto it = dstrb_parts_.begin(); it != dstrb_parts_.end(); ++it) {
		//cout << it->getName() << endl;
		for (std::size_t jj = 0; jj < it->getSize(); ++jj)
			duplicated_state_.SetBuffer( it->getDuplicatedStateIndex((int)jj),
					                     it->getData((int)jj) );
	}

	if (rank_ == numProcs_ - 1)
		for (auto it = shared_parts_.begin(); it != shared_parts_.end(); ++it) {
			//cout << it->getName() << endl;
			for (std::size_t jj = 0; jj < it->getSize(); ++jj)
				duplicated_state_.SetBuffer( it->getDuplicatedStateIndex((int)jj),
						                     it->getData((int)jj) );
		}

	duplicated_state_.Flush();

	if (rank_ == 0)
		std::cout << "[done]" << std::endl;

	return duplicated_state_;

}

/*
 *    This function essentially performs
 *    the opposite task of GetState.
 *    It takes the values from
 *    duplicated_state_ and updates the appropriate
 *    PHASTA Fortran arrays via the pointers
 *    stored in dstrb_parts_ and shared_parts_.
 */
void SimvascularVerdandiModel::StateUpdated() {

	if (rank_ == 0)
		std::cout << "setting state ";

	for(auto it = dstrb_parts_.begin(); it != dstrb_parts_.end(); ++it)
		for (std::size_t jj = 0; jj < it->getSize(); ++jj)
			it->setData( (int)jj, duplicated_state_( it->getDuplicatedStateIndex((int)jj) ) );

	double* tempArray = new double[shared_parts_size_];
	int tcounter = 0;

	if (rank_ == numProcs_ - 1) {
		for (auto it = shared_parts_.begin(); it != shared_parts_.end(); ++it)
			for (std::size_t jj = 0; jj < it->getSize(); jj++) {
				it->setData( (int)jj, duplicated_state_( it->getDuplicatedStateIndex((int)jj) ) );

				tempArray[tcounter++] = *it->getDataPointer((int)jj);
			}

		// communicate values from the last processor to the other processors
		MPI_Bcast(tempArray, shared_parts_size_, MPI_DOUBLE, numProcs_-1, iNewComm_C_);
	}
	else {

		MPI_Bcast(tempArray, shared_parts_size_, MPI_DOUBLE, numProcs_-1, iNewComm_C_);

		tcounter = 0;
		for (auto it = shared_parts_.begin(); it != shared_parts_.end(); ++it)
			for (std::size_t jj = 0; jj < it->getSize(); jj++)
				*it->getDataPointer((int)jj) = tempArray[tcounter++];
	}

	delete [] tempArray;

	// since we have set the state only on the master image nodes, we need to
	// copy state values from the master image to the interprocessor boundary nodes
	if (numProcs_ > 1)
		estim_helpers_setstate_comm();

    if (rank_ == 0)
		std::cout << "[done]" << std::endl;

}

////////////
// ERRORS //
////////////

/*! Returns a decomposition of the initial state error covariance matrix as a product \f$LUL^T\f$.
 */
/*!
      \param[out] L the matrix \f$L\f$.
      \param[out] U the matrix \f$U\f$.
 */
template <class L_matrix, class U_matrix>
void SimvascularVerdandiModel::GetStateErrorVarianceSqrt(L_matrix& L, U_matrix& U) {

	// L is distributed across processes (since one of its dimensions is Nstate)
	// but U is on every processes (dimensions of Nparameter)

	L.Reallocate(Nstate_, Nreduced_, Nstate_local_);
	L.Zero(); // TODO

	int start_ind_local, end_ind_local;
	L.GetProcessorRowRange(start_ind_local, end_ind_local);

	if (rank_ == numProcs_ - 1) {// it's zero except in the last processor

		int ncounter = 0;

		for (auto it = shared_parts_.begin(); it != shared_parts_.end(); ++it)
			for (std::size_t jj = 0; jj < it->getSize(); jj++)
				if (it->getIsEstimated((int)jj)) {
					L.SetBuffer(start_ind_local + state_reduced_start_local_ + ncounter, ncounter, double(1));
					ncounter++;
				}

	}
	L.Flush();
	//L.Print();

	U.Reallocate(Nreduced_, Nreduced_);
	U.Zero();

	for (int i = 0; i < Nreduced_; i++)
		U(i, i) = double(double(1) / state_error_variance_value_[i]);

	//U.Print();

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
	{
		stringstream s_temp;
		s_temp << this->rank_;
		string outname = "U_initial-"+s_temp.str()+".dat";
		ofstream outfile(outname.c_str(), ofstream::app);
		outfile.precision(std::numeric_limits<double>::digits10);

		U.WriteText(outfile);
	}
#endif

	// Read error variance value from a configuration file
	// Mlt(1./error_variance_value, U_);
    //L.Fill(0.);
    //for (int i = 0; i < Nreduced_; i++)
    //    L(i + Nstate_ - Nreduced_, i) = 1.;

}

/*!
      \return The pointer to distributed state parts
 */
SimvascularAugStatePart& SimvascularVerdandiModel::GetAugStateDstrb(std::string name) const {

    return *dstrb_parts_map_.at(name);

}

/*!
      \return The pointer to shared state parts
 */
SimvascularAugStatePart& SimvascularVerdandiModel::GetAugStateShared(std::string name) const {

    return *shared_parts_map_.at(name);

}


/*
 *    Here, we write out to file
 *    only the part of augmented state that is estimated
 *    (i.e., the parameters).
 */
void SimvascularVerdandiModel::WriteEstimates() {
	// write down the parameter values in a file
	//int state_start, state_end;

	//duplicated_state_.GetProcessorRange(state_start, state_end);

	if (rank_ == numProcs_ -1) {

		int ncounter = 0;

		for (auto it = shared_parts_.begin(); it != shared_parts_.end(); ++it)
			for (std::size_t jj = 0; jj < it->getSize(); jj++)
				if (it->getIsEstimated((int)jj)) {
					Eoutfile_ << *it->getDataPointer((int)jj) << " ";
					ncounter++;
				}

		Eoutfile_ << std::endl;

	}
}


/*!
      \return The number of MPI processes.
 */
int SimvascularVerdandiModel::GetNumProcs() const {
    return numProcs_;
}


/*!
      \return The rank of the MPI process
 */
int SimvascularVerdandiModel::GetRank() const {
	return rank_;
}



/*!
      \return The name of this class
 */
string SimvascularVerdandiModel::GetName() const {

	return "SimvascularVerdandiModel";
}



/*
      \param[in] message the received message.
 */
void SimvascularVerdandiModel::Message(string message) {

	// Put here any processing you need.
}


}

#endif
