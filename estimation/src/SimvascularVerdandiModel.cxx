#ifndef VERDANDI_FILE_MODEL_SimvascularVerdandiModel_CXX

#include "SimvascularVerdandiModel.hxx"

//! Checks for existance of file
/*!
      \param[in] const std::string& name
 */
inline bool exists_test1 (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

//! Constructor
SimvascularAugStatePart::SimvascularAugStatePart()
:   c_(1.0)
{
	name_ = "undefined";
}

//! Destructor
SimvascularAugStatePart::~SimvascularAugStatePart() {

}

//! Initialize with name
/*!
      \param[in] const std::string& setname
 */
void SimvascularAugStatePart::Initialize(const std::string& setname) {
	name_ = setname;
}

//! Adds a pointer to the array of pointers
/*!
      \param[in] double* data_pointer
 */
void SimvascularAugStatePart::addDataPointer(double* data_pointer) {
	pointers_.push_back(data_pointer);
}

//! Adds to the array of boolean flags for the estimated variables
/*!
      \param[in] bool val
 */
void SimvascularAugStatePart::addIsEstimated(bool val) {
	is_estimated_.push_back(val);
}

//! Sets a specific pointer in the array of pointers
/*!
      \param[in] int position
      \param[in] double* data_pointer
 */
void SimvascularAugStatePart::setDataPointer(int position, double* data_pointer) {
	pointers_[position] = data_pointer;
}

//! Sets the data for a specific variable; re-parameterize if it is an estimated variable
/*!
      \param[in] int position
      \param[in] double val
 */
void SimvascularAugStatePart::setData(int position, double val) {
	is_estimated_[position] ? *(pointers_[position]) = pow(2.0,val) : *(pointers_[position]) = val;
}


//! Sets the estimated-variable-flag for a specific variable
/*!
      \param[in] int position
      \param[in] bool val
 */
void SimvascularAugStatePart::setIsEstimated(int position, bool val) {
	is_estimated_[position] = val;
}

//! Sets the constant that the output of getData is multiplied with
/*!
      \param[in] double val
 */
void SimvascularAugStatePart::setPremulConstant(double val) {
	c_ = val;
}

//! Returns the data pointer to a specific variable
/*!
      \param[in] int position
 */
double * SimvascularAugStatePart::getDataPointer(int position) {
	return pointers_[position];
}

//! Returns the value of a specific variable; re-parameterize if it is an estimated variable
/*!
      \param[in] int position
 */
double SimvascularAugStatePart::getData(int position) {
	return is_estimated_[position] ? log2(*(pointers_[position]))*c_ : (*pointers_[position])*c_;
}

//! Returns the estimated-variable-flag for a specific variable
/*!
      \param[in] int position
 */
bool SimvascularAugStatePart::getIsEstimated(int position) {
	return is_estimated_[position];
}

//! Returns the number of variables
std::size_t SimvascularAugStatePart::getSize() {
	return pointers_.size();
}

//! Returns the number of estimated variables
std::size_t SimvascularAugStatePart::getNumEstimated() {
	unsigned int sum_of_elems=0;
	for(std::vector<bool>::iterator j=is_estimated_.begin();j!=is_estimated_.end();++j)
	    sum_of_elems += *j;
	return sum_of_elems;
}

//! Returns the index of the first estimated variable
int SimvascularAugStatePart::getFirstEstimated() {
	std::vector<bool>::iterator j;
	int k;
	for(j=is_estimated_.begin(), k = 0;j!=is_estimated_.end();j++,k++) {
		if (*j) {
			return k;
		}
	}
	return -1;
}

//! Returns the name
std::string SimvascularAugStatePart::getName() {
    return name_;
}

//! Empty all the arrays
void SimvascularAugStatePart::Clear() {
	pointers_.clear();
	is_estimated_.clear();
	state_error_variance_value_.clear();
}

namespace Verdandi {


////////////////////////////////
// CONSTRUCTOR AND DESTRUCTOR //
////////////////////////////////


//! Constructor.
SimvascularVerdandiModel::SimvascularVerdandiModel()
:   gat(NULL),
 	Nreduced_(0),
 	Nstate_(0),
 	Nstate_local_(0),
 	state_reduced_start_local_(-1),
 	est_parts_size_(0),
 	nreduced_has_wall_parameters_(0),
 	nreduced_has_coupled_parameters_(0),
 	cp_rcr_estimate_resistance_(0),
 	cp_rcr_estimate_compliance_(0),
 	cp_rcr_estimate_prox_resistance_(0),
 	cp_rcr_estimate_pout_(0),
 	cp_rcr_estimate_pstates_(0),
 	rank_(0),
 	numProcs_(1),
 	state_error_variance_value_(1),
 	iNewComm_C_(MPI_COMM_WORLD)
{
}


//! Destructor.
SimvascularVerdandiModel::~SimvascularVerdandiModel() {
	Eoutfile_.close();
	//MPI_Finalize();
}


////////////////
// INITIALIZE //
////////////////


//! Initializes the model.
/*!
      \param[in] configuration_file configuration file.
 */
void SimvascularVerdandiModel::Initialize(string configuration_file) {

	VerdandiOps configuration(configuration_file);

	configuration.SetPrefix("simvascular_model.");
	configuration.Set("state_reduced_has_wall_parameters",nreduced_has_wall_parameters_);
	configuration.Set("state_reduced_has_coupled_parameters",nreduced_has_coupled_parameters_);

	if (nreduced_has_coupled_parameters_) {
		configuration.Set("RCR_parameters_info.estimate_resistance",cp_rcr_estimate_resistance_);
		configuration.Set("RCR_parameters_info.estimate_compliance",cp_rcr_estimate_compliance_);
		configuration.Set("RCR_parameters_info.estimate_prox_resistance",cp_rcr_estimate_prox_resistance_);
		configuration.Set("RCR_parameters_info.estimate_pout",cp_rcr_estimate_pout_);
		configuration.Set("RCR_parameters_info.estimate_pstates",cp_rcr_estimate_pstates_);

		configuration.Set("RCR_parameters_info.resistance_included",cp_rcr_include_resistance_);
		configuration.Set("RCR_parameters_info.compliance_included",cp_rcr_include_compliance_);
		configuration.Set("RCR_parameters_info.prox_resistance_included",cp_rcr_include_prox_resistance_);
	}

	configuration.Set("error_statistics.state_error_variance",state_error_variance_value_);

    Initialize();

}

//! Initializes the model.
/*!

 */
void SimvascularVerdandiModel::Initialize() {

	char pathToProcsCaseDir[100];

	// Get pointer to the single instance of PhGlobalArrayTransfer
	gat = PhGlobalArrayTransfer::Get();

	// save the communicator
	iNewComm_C_ = MPI_COMM_WORLD;
	newcom.iNewComm = MPI_Comm_c2f(iNewComm_C_); // modifies newcom in fortran common block

	MPI_Comm_size(iNewComm_C_, &numProcs_);
	MPI_Comm_rank(iNewComm_C_, &rank_);

	// read configuration file
	input_fform();

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

    // find nonzero entries of wall props
    for (int surfIdx = 0; surfIdx <= MAXSURF; surfIdx++)
    	if (nomodule.ValueListWallE[surfIdx] > 0)
    		WallEInd_.PushBack(surfIdx);

    // organize the arrays from PHASTA into the SimvascularAugStatePart class
    BuildAugmentedState();

	// Compute the local state size
	// perhaps this can be less hard-coded
	// NSD is 3
	// 4 comes from the 3 vel components and 1 pressure component

    /*
    Nreduced_ = 0;
	Nstate_local_ = (2 * 4) * conpar.nshguniq;

	if (nomodule.ideformwall > 0) {
		Nstate_local_ += NSD * conpar.nshguniq;

		if (nomodule.imeasdist > 0)
			Nstate_local_ += conpar.nshguniq * 2; // distance vector and pre-multiplied distance vector
	}

	if (rank_ == numProcs_ - 1) {

		// Add the lumped parameter states to the local state size of the last process
		//Nstate_local_ += grcrbccom.numGRCRSrfs * 2; // number of lumped parameter states

		Nstate_local_ += grcrbccom.numGRCRSrfs;

		if (!(nreduced_has_coupled_parameters_ && cp_rcr_estimate_pstates_))
			Nstate_local_ += grcrbccom.numGRCRSrfs;

		state_reduced_start_local_ = Nstate_local_; // index for first element of reduced state

		if (nreduced_has_coupled_parameters_ && cp_rcr_estimate_pstates_)
			Nstate_local_ += grcrbccom.numGRCRSrfs;

	}

	// Compute the reduced state size
	// (the number of variables to be estimated)
    if (nreduced_has_wall_parameters_  && nomodule.ideformwall > 0) {
    	Nreduced_ += nomodule.numWallRegions;

    	// Add the number of wall parameters to the local state size of the last process
    	if (rank_ == numProcs_ - 1)
    		Nstate_local_ += nomodule.numWallRegions;
    }

    if (nreduced_has_coupled_parameters_) {

    	if (cp_rcr_estimate_pstates_)
    		Nreduced_ += grcrbccom.numGRCRSrfs; // account for reduced-order P state

    	if (cp_rcr_estimate_compliance_) {

    		for (unsigned int kk = 0; kk < cp_rcr_include_compliance_.size(); kk++) {
    			if (cp_rcr_include_compliance_[kk]) {
    				Nreduced_++;
    				if (rank_ == numProcs_ - 1)
    					Nstate_local_++;
    			}
    		}
    	}

    	if (cp_rcr_estimate_resistance_) {

    		for (unsigned int kk = 0; kk < cp_rcr_include_resistance_.size(); kk++) {
    			if (cp_rcr_include_resistance_[kk]) {
    				Nreduced_++;
    				if (rank_ == numProcs_ - 1)
    					Nstate_local_++;
    			}
            }
    	}

    	if (cp_rcr_estimate_prox_resistance_) {

    		for (unsigned int kk = 0; kk < cp_rcr_include_prox_resistance_.size(); kk++) {
    			if (cp_rcr_include_prox_resistance_[kk]) {
    				Nreduced_++;
    				if (rank_ == numProcs_ - 1)
    					Nstate_local_++;
    			}
    		}
    	}

    	if (cp_rcr_estimate_pout_) {
    		Nreduced_++;
    		if (rank_ == numProcs_ - 1)
    			Nstate_local_++;
    	}

    }
    */

    int estimated_total = 0, state_total = 0;

    std::cout << "at rank " << rank_ << " number of states " << state_parts_.size() << std::endl;

    for (std::size_t kk = 0; kk < state_parts_.size(); kk++) {

    	std::cout << "at rank " << rank_ << " " << state_parts_[(int)kk].getName().c_str() << " " << state_parts_[(int)kk].getSize() << std::endl;

    	state_total += (int)state_parts_[(int)kk].getSize();

    }

    if (rank_ == numProcs_ - 1) {

    	for (std::size_t kk = 0; kk < est_parts_.size(); kk++) {

    		std::cout << "at rank " << rank_ << " " << est_parts_[(int)kk].getName().c_str() << " " << est_parts_[(int)kk].getSize() << std::endl;

    		int first_estimated = est_parts_[(int)kk].getFirstEstimated();

    		if (first_estimated >= 0 && state_reduced_start_local_ < 0)
    			state_reduced_start_local_ = (int)state_total+first_estimated;

    		estimated_total += (int)est_parts_[(int)kk].getNumEstimated();

    		state_total += (int)est_parts_[(int)kk].getSize();

    	}
    }

    for (std::size_t kk = 0; kk < est_parts_.size(); kk++)
    	est_parts_size_ += est_parts_[(int)kk].getSize();

    Nstate_local_ = state_total;

    Nreduced_ = estimated_total;

    // we want Nreduced_ on every proc
    MPI_Allreduce(&Nreduced_, &Nreduced_, 1, MPI_INT, MPI_SUM, iNewComm_C_);

    std::cout << "at rank " << rank_ << " # of distributed states: " << Nstate_local_ << std::endl;
    std::cout << "at rank " << rank_ << " # of shared states: " << est_parts_size_ << std::endl;
    std::cout << "at rank " << rank_ << " # of estimated states: " << Nreduced_ << std::endl;
    std::cout << "at rank " << rank_ << " state_reduced_start " << state_reduced_start_local_ << std::endl;

	// Compute the global state size
	MPI_Allreduce(&Nstate_local_, &Nstate_, 1, MPI_INT, MPI_SUM, iNewComm_C_);

	// Make sure the state array is distributed with the correct number
	// of elements per processor that is given by Nstate_local_
	duplicated_state_.Reallocate(Nstate_, Nstate_local_);

    // file output for initial step if starting from scratch
    std::string Efilename = "estimated_reduced_state.dat";

    if (exists_test1(Efilename)) { // still need to finish this
        // The file exists, and is open for input
    	std::cout << "we can start from the previous estimates" << std::endl;
    } else {
    	std::cout << "we need to write the new file" << std::endl;
    }

    if (rank_ == numProcs_ - 1) {

    	Eoutfile_.open(Efilename.c_str(),ofstream::app);
    	Eoutfile_ << std::scientific << std::setprecision( std::numeric_limits<double>::digits10 );

    }

}


//! Organizes the array pointers from simvascular
void SimvascularVerdandiModel::BuildAugmentedState() {

	SimvascularAugStatePart state_part;

	// velocity and pressure field
	state_part.Initialize("solution (v1,v2,v3,p)");
	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {
		int actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];
		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond
			state_part.addDataPointer(&(gat->global_yold_ptr[varIdx * conpar.nshg + actualIdx-1]));
			state_part.addIsEstimated(0);
		}
	}
	state_parts_.push_back(state_part);
	state_part.Clear();

	// acceleration fields
	state_part.Initialize("time derivative of solution (v1',v2',v3',p')");
	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {
		int actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];
		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond
			state_part.addDataPointer(&(gat->global_acold_ptr[varIdx * conpar.nshg + actualIdx-1]));
			state_part.addIsEstimated(0);
		}
	}
	state_parts_.push_back(state_part);
	state_part.Clear();

	if (nomodule.ideformwall > 0) {

		// displacement field
		state_part.Initialize("displacement (u1,u2,u3)");
		for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {
			int actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];
			for(int varIdx=0; varIdx < 3; varIdx++) {
				state_part.addDataPointer(&(gat->global_uold_ptr[varIdx * conpar.nshg + actualIdx-1]));
				state_part.addIsEstimated(0);
			}
		}
		state_parts_.push_back(state_part);
		state_part.Clear();

		if (nomodule.imeasdist > 0) {

			// these distance fields are stored in the state strictly for convenience
			// distance field
			state_part.Initialize("distance (d)");
			for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {
				int actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];
				for(int varIdx=0; varIdx < 1; varIdx++) {
					state_part.addDataPointer(&(gat->global_xdist_ptr[varIdx * conpar.nshg + actualIdx-1]));
					state_part.addIsEstimated(0);
					state_part.setPremulConstant(-1.0);
				}
			}
			state_parts_.push_back(state_part);
			state_part.Clear();

			// mass matrix multiplied by distance
			state_part.Initialize("M * distance (M*d)");
			for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {
				int actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];
				for(int varIdx=0; varIdx < 1; varIdx++) {
					state_part.addDataPointer(&(gat->global_df_fem_ptr[varIdx * conpar.nshg + actualIdx-1]));
					state_part.addIsEstimated(0);
					state_part.setPremulConstant(-1.0);
				}
			}
			state_parts_.push_back(state_part);
            state_part.Clear();

		}

	}


	// this is where the estimated parts of the augmented state will go

	// lumped parameter state (for RCR, just P_current)
	if (grcrbccom.numGRCRSrfs > 0) {
		state_part.Initialize("RCR P_current");
		for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++) {
			state_part.addDataPointer(&(gat->global_lumped_parameter_P[surfIdx]));
			if (cp_rcr_estimate_pstates_)
				state_part.addIsEstimated(1);
			else
				state_part.addIsEstimated(0);
		}
		est_parts_.push_back(state_part);
		state_part.Clear();
	}

	// for vessel wall regional properties (E)
	if (nreduced_has_wall_parameters_ && nomodule.ideformwall > 0) {
		state_part.Initialize("E");
		for(int parIdx = 0; parIdx < nomodule.numWallRegions; parIdx++) {
			state_part.addDataPointer(&(nomodule.ValueListWallE[WallEInd_(parIdx)])); //note not from global array transfer (hmm...)
			state_part.addIsEstimated(1);
		}
		est_parts_.push_back(state_part);
		state_part.Clear();
	}

	// for lumped parameter model parameters
	// for RCR
	if (nreduced_has_coupled_parameters_ && grcrbccom.numGRCRSrfs > 0) {

		// for RCR C
		if (cp_rcr_estimate_compliance_) {

			state_part.Initialize("RCR C");
			for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++) {

				if (cp_rcr_include_compliance_[parIdx]) {
					state_part.addDataPointer(&(gat->global_lumped_parameter_params[parIdx*3+1]));
					state_part.addIsEstimated(1);
				}

			}
			est_parts_.push_back(state_part);
			state_part.Clear();

		}

		// for RCR R2
		if (cp_rcr_estimate_resistance_) {

			state_part.Initialize("RCR R2");
			for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++) {

				if (cp_rcr_include_resistance_[parIdx]) {
					state_part.addDataPointer(&(gat->global_lumped_parameter_params[parIdx*3+2]));
					state_part.addIsEstimated(1);
				}

			}
			est_parts_.push_back(state_part);
			state_part.Clear();

		}

		// for RCR R1
		if (cp_rcr_estimate_prox_resistance_) {

			state_part.Initialize("RCR R1");
			for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++) {

				if (cp_rcr_include_prox_resistance_[parIdx]) {
					state_part.addDataPointer(&(gat->global_lumped_parameter_params[parIdx*3+0]));
					state_part.addIsEstimated(1);
				}

			}
			est_parts_.push_back(state_part);
			state_part.Clear();

		}

		// for RCR Pdist
		if (cp_rcr_estimate_pout_) {

			state_part.Initialize("RCR Pdist");

			state_part.addDataPointer(&(gat->global_lumped_parameter_pout[0]));
			state_part.addIsEstimated(1);

			est_parts_.push_back(state_part);
			state_part.Clear();

		}

	}

}

//! Initializes the first time step for the model.
void SimvascularVerdandiModel::InitializeFirstStep() {

}

//! Initializes the current time step for the model.
void SimvascularVerdandiModel::InitializeStep() {

}

//! Finalizes the current time step for the model.
void SimvascularVerdandiModel::Finalize() {

	itrdrv_finalize();

}


////////////////
// PROCESSING //
////////////////


//! Advances one step forward in time.
/*! \f[x^f_{h+1} = \mathcal{M}_h(x^a_h, p_h)\,.\f] */
void SimvascularVerdandiModel::Forward() {

	itrdrv_iter_init();

	itrdrv_iter_step();

	itrdrv_iter_finalize();

}

//! Finalizes the current time step
//! meant to be called after the mean state with innovation is set
void SimvascularVerdandiModel::FinalizeStep() {

	itrdrv_iter_finalize(); // routines that allow moving to the next step

}



//! Checks whether the model has finished.
/*!
      \return True if the simulation is done, false otherwise.
 */
bool SimvascularVerdandiModel::HasFinished() const {

	return timdat.istep >= inpdat.nstep[0];
}


///////////////
// OPERATORS //
///////////////


//! Applies the model to a given vector.
/*! The current state of the model is modified.
      \param[in] x a vector.
      \param[in] forward Boolean to indicate if the model has to go on to the
      next step.
      \param[in] preserve_state Boolean to indicate if the model state has to
      be preserved.
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


//! Returns the current time.
/*!
      \return The current time.
 */
double SimvascularVerdandiModel::GetTime() const {

	return (double)(timdat.lstep);
}


//! Sets the time of the model to a given time.
/*!
      \param[in] time a given time.
 */
void SimvascularVerdandiModel::SetTime(double time) {
	throw ErrorUndefined("SimvascularVerdandiModel"
			"::SetTime(double time)");
}


//! Returns the state vector size.
/*!
      \return The state vector size.
 */
int SimvascularVerdandiModel::GetNstate() const {
	return Nstate_;
}


//! Returns the size of the local (on-processor) state vector.
/*!
      \return The size of the local state vector.
 */
int SimvascularVerdandiModel::GetLocalNstate() const {
	return Nstate_local_; // this is used in reallocate routine for petsc matrices in ROUKF
}

//! Returns the starting index of the local reduced state vector
/*!
      \return The starting index of the local reduced state vector.
 */
int SimvascularVerdandiModel::GetLocalReducedStart() const {
	return state_reduced_start_local_;
}


//! copies internal state to "duplicated state" and returns a reference to it
/*!
      \return: a reference to "duplicated state" vector
 */
SimvascularVerdandiModel::state& SimvascularVerdandiModel::GetState() {

	int err;
	double val;
	int state_start, state_end, icounter;

	int actualIdx;

	if (rank_ == 0)
		std::cout << "getting state ";

	duplicated_state_.GetProcessorRange(state_start, state_end);

	icounter = state_start;

	for (std::size_t kk = 0; kk < state_parts_.size(); kk++) {
		for (std::size_t jj = 0; jj < state_parts_[(int)kk].getSize(); jj++) {
			duplicated_state_.SetBuffer(icounter++,state_parts_[(int)kk].getData((int)jj));
		}
	}

	if (rank_ == numProcs_ - 1) {
		for (std::size_t kk = 0; kk < est_parts_.size(); kk++) {
			for (std::size_t jj = 0; jj < est_parts_[(int)kk].getSize(); jj++) {
				duplicated_state_.SetBuffer(icounter++,est_parts_[(int)kk].getData((int)jj));
			}
		}
	}

	/*
	// Get the "solution field"
	// The state elements must be distributed such that there are no duplicate entries across processors
	// The list of unique nodes (master image) is subset of the nodes that are on the local processor

	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

		actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond

			duplicated_state_.SetBuffer(icounter++,(gat->global_yold_ptr)[varIdx * conpar.nshg + actualIdx-1]);

		}

	}

	// get the "acceleration field"
	// in the flowsolver this is acold

	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

		actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond

			duplicated_state_.SetBuffer(icounter++,(gat->global_acold_ptr)[varIdx * conpar.nshg + actualIdx-1]);

		}

	}

	// get the "displacement field"
	// in the flowsolver this is uold

	if (nomodule.ideformwall > 0) {

		for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

			actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

			for(int varIdx=0; varIdx < 3; varIdx++) { // 3 dofs

				duplicated_state_.SetBuffer(icounter++,(gat->global_uold_ptr)[varIdx * conpar.nshg + actualIdx-1]);

			}

		}

	}

	// get the "distance field"
	// even though the distances are not true state variables
	// (i.e. they are not used to in the current time step to update the
	// state variables for the next time step), it is convenient to store the distances
	// in the state vector that is passed to the roukf filter driver so that they may be accessed
	// by the observation manager (for the distance field innovation)

	if (nomodule.ideformwall > 0 && nomodule.imeasdist > 0) {

		for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

			actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

			// in the observation manager, the innovation takes the general form of G(X) = Z-H(X)
			// for practical purpooses, the applyoperator function in
			// the observation manager assigns H(X) to the signed distance and zero to Z;

			duplicated_state_.SetBuffer(icounter++,-(gat->global_xdist_ptr)[actualIdx-1]);

		}

		// get the distance field premultiplied by a finite element 'mass' matrix
		// on the wall boundary
		// i.e. the effect of the mass matrix part of the measurement covariance is
		// already included here

		for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

			actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

			duplicated_state_.SetBuffer(icounter++,-(gat->global_df_fem_ptr)[actualIdx-1]);
			//duplicated_state_.SetBuffer(icounter++,-(gat->global_xdist_ptr)[actualIdx-1]);

		}

	}


	// get the P states for the lumped parameter model
	// note that this part of the state is only on the last processor

	if (rank_ == numProcs_ - 1) {

		//for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++)
		//	duplicated_state_.SetBuffer(icounter++,(gat->global_lumped_parameter_Q)[surfIdx]);

		for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++)
			duplicated_state_.SetBuffer(icounter++,log2(gat->global_lumped_parameter_P[surfIdx]));

	}

	// get the parameter values
	// note that this part of the state is only on the last processor

	if (rank_ == numProcs_ - 1) {

		if (nreduced_has_wall_parameters_ && nomodule.ideformwall > 0) {

			//val = nomodule.evw;
			//duplicated_state_.SetBuffer(icounter++,log2(val));

			for(int parIdx = 0; parIdx < nomodule.numWallRegions; parIdx++) {

				val = nomodule.ValueListWallE[WallEInd_(parIdx)];

				duplicated_state_.SetBuffer(icounter++,log2(val));

			}
			//std::cout << "[get] rank: " << rank_ << " val: " << log2(val) << std::endl;
		}

		if (nreduced_has_coupled_parameters_ && grcrbccom.numGRCRSrfs > 0) {

			if (cp_rcr_estimate_compliance_)
				for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++) {

					if (cp_rcr_include_compliance_[parIdx]) {
						val = gat->global_lumped_parameter_params[parIdx*3+1];

						duplicated_state_.SetBuffer(icounter++,log2(val));

						//				std::cout << "outlet " << parIdx << " ";
						//				std::cout << gat->global_lumped_parameter_params[parIdx*3+0] << " ";
						//			    std::cout << gat->global_lumped_parameter_params[parIdx*3+1] << " ";
						//				std::cout << gat->global_lumped_parameter_params[parIdx*3+2] << " " << std::endl;
					}

				}

			if (cp_rcr_estimate_resistance_)
				for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++) {

					if (cp_rcr_include_resistance_[parIdx]) {

						val = gat->global_lumped_parameter_params[parIdx*3+2];

						duplicated_state_.SetBuffer(icounter++,log2(val));

					}

				}

			if (cp_rcr_estimate_prox_resistance_)
				for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++) {

					if (cp_rcr_include_prox_resistance_[parIdx]) {

						val = gat->global_lumped_parameter_params[parIdx*3+0];

						duplicated_state_.SetBuffer(icounter++,log2(val));

					}

				}

			if (cp_rcr_estimate_pout_) {
				val = gat->global_lumped_parameter_pout[0];
				duplicated_state_.SetBuffer(icounter++,log2(val));
			}

		}

	}*/

	duplicated_state_.Flush();

	if (rank_ == 0)
		std::cout << "[done]" << std::endl;

	return duplicated_state_;

}

//! Updates the internal state from the duplicated state
void SimvascularVerdandiModel::StateUpdated() {

	int err;
	double val;
	int state_start, state_end, icounter;

	int actualIdx;

	if (rank_ == 0)
		std::cout << "setting state ";

	duplicated_state_.GetProcessorRange(state_start, state_end);

	icounter = state_start;


	for (std::size_t kk = 0; kk < state_parts_.size(); kk++) {
		for (std::size_t jj = 0; jj < state_parts_[(int)kk].getSize(); jj++) {
			state_parts_[(int)kk].setData((int)jj,duplicated_state_(icounter++));
		}
	}

	double* tempArray = new double[est_parts_size_];
	int tcounter = 0;

	if (rank_ == numProcs_ - 1) {
		for (std::size_t kk = 0; kk < est_parts_.size(); kk++) {
			for (std::size_t jj = 0; jj < est_parts_[(int)kk].getSize(); jj++) {
				est_parts_[(int)kk].setData((int)jj,duplicated_state_(icounter++));
				tempArray[tcounter++] = *est_parts_[(int)kk].getDataPointer((int)jj);
			}
		}

		MPI_Bcast(tempArray, est_parts_size_, MPI_DOUBLE, numProcs_-1, iNewComm_C_);
	}
	else {

		MPI_Bcast(tempArray, est_parts_size_, MPI_DOUBLE, numProcs_-1, iNewComm_C_);

		tcounter = 0;
		for (std::size_t kk = 0; kk < est_parts_.size(); kk++)
			for (std::size_t jj = 0; jj < est_parts_[(int)kk].getSize(); jj++)
				*est_parts_[(int)kk].getDataPointer((int)jj) = tempArray[tcounter++];
	}

	delete [] tempArray;

    /*
	// set the "solution field"
	// The list of unique nodes (master image) is subset of the nodes that are on the local processor

	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

		actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond

			(gat->global_yold_ptr)[varIdx * conpar.nshg + actualIdx-1] = duplicated_state_(icounter++);

		}

	}

	// set the "acceleration field"
	// in the flowsolver this is acold

	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

		actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond

			(gat->global_acold_ptr)[varIdx * conpar.nshg + actualIdx-1] = duplicated_state_(icounter++);

		}

	}

	// set the "displacement field"
	// in the flowsolver this is uold

	if (nomodule.ideformwall > 0) {

		for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

			actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

			for(int varIdx=0; varIdx < 3; varIdx++) { // 3 dofs

				(gat->global_uold_ptr)[varIdx * conpar.nshg + actualIdx-1] = duplicated_state_(icounter++);

			}

		}

	}

	// increment the counter to skip the "distance field"
	// since the distance field does not get used to update the state
	// we don't need to update it internally, but the counter
	// needs to be incremented

	if (nomodule.imeasdist > 0) {

		for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {
			actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

			//icounter++;
			(gat->global_xdist_ptr)[actualIdx-1] = duplicated_state_(icounter++);
			//icounter++; // skips the space occupied by the second finite element distance vector
		}

		for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {
			actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

			(gat->global_df_fem_ptr)[actualIdx-1] = duplicated_state_(icounter++);
		}

	}*/

	/*

	// set the P,Q states for the lumped parameter model

	if (rank_ == numProcs_ - 1) {

		// setting the Q state isn't strictly necessary
		// as Q is recomputed from the velocity field
		//for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++) {

		//	(gat->global_lumped_parameter_Q)[surfIdx] = duplicated_state_(icounter++);

		//}

		for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++) {

			//(gat->global_lumped_parameter_P)[surfIdx] = pow(2.0,duplicated_state_(icounter++));
			(gat->global_lumped_parameter_P)[surfIdx] = duplicated_state_(icounter++);

		}
	}

	MPI_Bcast(gat->global_lumped_parameter_P,grcrbccom.numGRCRSrfs, MPI_DOUBLE, numProcs_-1, iNewComm_C_);
	//MPI_Bcast(gat->global_lumped_parameter_Q,grcrbccom.numGRCRSrfs, MPI_DOUBLE, numProcs_-1, iNewComm_C_);
	//std::cout << "average distal P at rank " << rank_ << " : ";

//	for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++) {
//		std::cout << (gat->global_lumped_parameter_P)[surfIdx] << " ";
//	}
//	std::cout << std::endl;


	// set the values of the reduced state
	// note that this part of the state is only on the last processor
	// so we need to broadcast the parameter value to every processor

	if (nreduced_has_wall_parameters_ && nomodule.ideformwall > 0) {

		double* tempArray = new double[nomodule.numWallRegions];

		if (rank_ == numProcs_ - 1)
			for(int parIdx = 0; parIdx < nomodule.numWallRegions; parIdx++)
				tempArray[parIdx] = pow(2.0,duplicated_state_(icounter++));

		//MPI_Bcast(&val, 1, MPI_DOUBLE, numProcs_ - 1, iNewComm_C_);
		//nomodule.evw = val;

		MPI_Bcast(tempArray, nomodule.numWallRegions, MPI_DOUBLE, numProcs_ - 1, iNewComm_C_);

		for(int parIdx = 0; parIdx < nomodule.numWallRegions; parIdx++)
			nomodule.ValueListWallE[WallEInd_(parIdx)] = tempArray[parIdx];

		delete [] tempArray;

	}

	if (nreduced_has_coupled_parameters_ && grcrbccom.numGRCRSrfs > 0) {

		double* tempArray = new double[grcrbccom.numGRCRSrfs];

		if (cp_rcr_estimate_compliance_) {
			if (rank_ == numProcs_ - 1)
				for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++)
					if (cp_rcr_include_compliance_[parIdx]) {
						tempArray[parIdx] = pow(2.0,duplicated_state_(icounter++));
					} else {
						tempArray[parIdx] = 0.0;
					}

			MPI_Bcast(tempArray, grcrbccom.numGRCRSrfs, MPI_DOUBLE, numProcs_ - 1, iNewComm_C_);

			for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++)
				if (cp_rcr_include_compliance_[parIdx])
					gat->global_lumped_parameter_params[parIdx*3+1] = tempArray[parIdx];
		}

		if (cp_rcr_estimate_resistance_) {
			if (rank_ == numProcs_ - 1)
				for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++)
					if (cp_rcr_include_resistance_[parIdx]) {
						tempArray[parIdx] = pow(2.0,duplicated_state_(icounter++));
					} else {
						tempArray[parIdx] = 0.0;
					}

			MPI_Bcast(tempArray, grcrbccom.numGRCRSrfs, MPI_DOUBLE, numProcs_ - 1, iNewComm_C_);

			for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++)
				if (cp_rcr_include_resistance_[parIdx])
					gat->global_lumped_parameter_params[parIdx*3+2] = tempArray[parIdx];
		}

		if (cp_rcr_estimate_prox_resistance_) {
			if (rank_ == numProcs_ - 1)
				for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++)
					if (cp_rcr_include_prox_resistance_[parIdx]) {
						tempArray[parIdx] = pow(2.0,duplicated_state_(icounter++));
					} else {
						tempArray[parIdx] = 0.0;
					}

			MPI_Bcast(tempArray, grcrbccom.numGRCRSrfs, MPI_DOUBLE, numProcs_ - 1, iNewComm_C_);

			for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++)
				if (cp_rcr_include_prox_resistance_[parIdx])
					gat->global_lumped_parameter_params[parIdx*3+0] = tempArray[parIdx];
		}

		if (cp_rcr_estimate_pout_) {
			if (rank_ == numProcs_ - 1) {
				tempArray[0] = pow(2.0,duplicated_state_(icounter++));
			}

			MPI_Bcast(tempArray, 1, MPI_DOUBLE, numProcs_ - 1, iNewComm_C_);

			gat->global_lumped_parameter_pout[0] = tempArray[0];
		}

		delete [] tempArray;
	}
    */

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

/*! Returns a decomposition of the initial state error covariance matrix (\f$B\f$)
      as a product \f$LUL^T\f$.
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

		for (std::size_t kk = 0; kk < est_parts_.size(); kk++) {
			for (std::size_t jj = 0; jj < est_parts_[(int)kk].getSize(); jj++) {
				if (est_parts_[(int)kk].getIsEstimated((int)jj)) {
					L.SetBuffer(start_ind_local + state_reduced_start_local_ + ncounter, ncounter, double(1));
					ncounter++;
				}
			}
		}

		/*
		if (nreduced_has_coupled_parameters_ && cp_rcr_estimate_pstates_) { // for the reduced-order P-states
			for (int i = 0; i < grcrbccom.numGRCRSrfs; i++) {
				L.SetBuffer(start_ind_local + state_reduced_start_local_ + ncounter, ncounter, double(1));
				ncounter++;
			}
		}

		if (nreduced_has_wall_parameters_)
			for (int i = 0; i < nomodule.numWallRegions; i++) {
				L.SetBuffer(start_ind_local + state_reduced_start_local_ + ncounter, ncounter, double(1));
				ncounter++;
			}

		if (nreduced_has_coupled_parameters_) {

			if (cp_rcr_estimate_compliance_)
				for (int i = 0; i < grcrbccom.numGRCRSrfs; i++) {
					if (cp_rcr_include_compliance_[i]) {
						L.SetBuffer(start_ind_local + state_reduced_start_local_ + ncounter, ncounter, double(1));
						ncounter++;
					}
				}

			if (cp_rcr_estimate_resistance_)
				for (int i = 0; i < grcrbccom.numGRCRSrfs; i++) {
					if (cp_rcr_include_resistance_[i]) {
						L.SetBuffer(start_ind_local + state_reduced_start_local_ + ncounter, ncounter, double(1));
						ncounter++;
					}
				}

			if (cp_rcr_estimate_prox_resistance_)
				for (int i = 0; i < grcrbccom.numGRCRSrfs; i++) {
					if (cp_rcr_include_prox_resistance_[i]) {
						L.SetBuffer(start_ind_local + state_reduced_start_local_ + ncounter, ncounter, double(1));
						ncounter++;
					}
				}

			if (cp_rcr_estimate_pout_){
				L.SetBuffer(start_ind_local + state_reduced_start_local_ + ncounter, ncounter, double(1));
				ncounter++;
			}

		}
		*/

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

//! Write only the estimated parameters to file
void SimvascularVerdandiModel::WriteEstimates() {
	// write down the parameter values in a file
	int state_start, state_end;

	duplicated_state_.GetProcessorRange(state_start, state_end);

	if (rank_ == numProcs_ -1) {

		int ncounter = 0;

		for (std::size_t kk = 0; kk < est_parts_.size(); kk++) {
			for (std::size_t jj = 0; jj < est_parts_[(int)kk].getSize(); jj++) {
				if (est_parts_[(int)kk].getIsEstimated((int)jj)) {
					Eoutfile_ << *est_parts_[(int)kk].getDataPointer((int)jj) << " ";
					ncounter++;
				}
			}
		}

		/*
		// reduced-order P-states
		if (nreduced_has_coupled_parameters_ && cp_rcr_estimate_pstates_ && grcrbccom.numGRCRSrfs > 0) {
			for (int kk = 0; kk < grcrbccom.numGRCRSrfs; kk++) {
				Eoutfile_ << pow(2.0,duplicated_state_(state_start + state_reduced_start_local_ + ncounter) ) << " ";
				ncounter++;
			}
		}

		// wall parameters
		if (nreduced_has_wall_parameters_ && nomodule.ideformwall > 0) {
			for (int kk = 0; kk < nomodule.numWallRegions; kk++) {
				Eoutfile_ << pow(2.0,duplicated_state_(state_start + state_reduced_start_local_ + ncounter) ) << " ";
				ncounter++;
			}
		}

		// reduced-order parameters
		if (nreduced_has_coupled_parameters_ && grcrbccom.numGRCRSrfs > 0) {

			if (cp_rcr_estimate_compliance_)
				for (int kk = 0; kk < grcrbccom.numGRCRSrfs; kk++) {
					if (cp_rcr_include_compliance_[kk]) {
						Eoutfile_ << pow(2.0,duplicated_state_(state_start + state_reduced_start_local_ + ncounter) ) << " ";
						ncounter++;
					}
				}

			if (cp_rcr_estimate_resistance_)
				for (int kk = 0; kk < grcrbccom.numGRCRSrfs; kk++) {
					if (cp_rcr_include_resistance_[kk]) {
						Eoutfile_ << pow(2.0,duplicated_state_(state_start + state_reduced_start_local_ + ncounter) ) << " ";
						ncounter++;
					}
				}

			if (cp_rcr_estimate_prox_resistance_)
				for (int kk = 0; kk < grcrbccom.numGRCRSrfs; kk++) {
					if (cp_rcr_include_prox_resistance_[kk]) {
						Eoutfile_ << pow(2.0,duplicated_state_(state_start + state_reduced_start_local_ + ncounter) ) << " ";
						ncounter++;
					}
				}

			if (cp_rcr_estimate_pout_) {
				Eoutfile_ << pow(2.0,duplicated_state_(state_start + state_reduced_start_local_ + ncounter) ) << " ";
				ncounter++;
			}

		}
		*/

		Eoutfile_ << std::endl;
	}
}

//! Returns the number of MPI
/*!
      \return The name of the class.
 */
int SimvascularVerdandiModel::GetNumProcs() const {
    return numProcs_;
}

//! Returns the name of the class.
/*!
      \return The name of the class.
 */
int SimvascularVerdandiModel::GetRank() const {
	return rank_;
}


//! Returns the name of the class.
/*!
      \return The name of the class.
 */
string SimvascularVerdandiModel::GetName() const {

	return "SimvascularVerdandiModel";
}


//! Receives and handles a message.
/*
      \param[in] message the received message.
 */
void SimvascularVerdandiModel::Message(string message) {

	// Put here any processing you need.
}


}

#define VERDANDI_FILE_MODEL_SimvascularVerdandiModel_CXX
#endif
