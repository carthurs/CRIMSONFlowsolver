#ifndef VERDANDI_FILE_MODEL_SimvascularVerdandiModel_CXX

#include "SimvascularVerdandiModel.hxx"

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
 	nreduced_has_wall_parameters_(0),
 	nreduced_has_coupled_parameters_(0),
 	rank_(0),
 	numProcs_(1),
 	state_error_variance_value_(1),
 	iNewComm_C_(MPI_COMM_WORLD)
{
}


//! Destructor.
SimvascularVerdandiModel::~SimvascularVerdandiModel() {
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
	configuration.Set("state_reduced_has_wall_parameters",
			nreduced_has_wall_parameters_);
	configuration.Set("state_reduced_has_coupled_parameters",
			nreduced_has_coupled_parameters_);


	if (nreduced_has_coupled_parameters_) {
		configuration.Set("RCR_parameters_info.estimate_resistance",
				cp_rcr_estimate_resistance_);
		configuration.Set("RCR_parameters_info.estimate_compliance",
				cp_rcr_estimate_compliance_);

		configuration.Set("RCR_parameters_info.resistance_included",
				cp_rcr_include_resistance_);
		configuration.Set("RCR_parameters_info.compliance_included",
				cp_rcr_include_compliance_);

//		configuration.Set("RCR_parameters_info.face_grouping",
//						cp_rcr_face_grouping_);
	}

	configuration.Set("error_statistics.state_error_variance",
			state_error_variance_value_);

	//cout << "assimilation file: " << configuration.GetFilePath() << endl;

    Initialize();

}

//! Initializes the model.
/*!

 */
void SimvascularVerdandiModel::Initialize() {

	char pathToProcsCaseDir[100];

//	MPI_Comm newcomm;
//	int color,key;
//	int numparticles = 1;
//	int numprocs_perparticle;

	//MPI_Init(&argc,&argv);

	// create MPI communicator for one simulation
	// by splitting MPI_COMM_WORLD
	// this functionality isn't used for now as
	// the "numparticles" always is 1
//	numprocs_perparticle = numProcsTotal / numparticles;
//	color = rank_ / numprocs_perparticle;
//	key = 0;
//	MPI_Comm_split(MPI_COMM_WORLD,color,key,&newcomm);

	// Initialize phSolver
	//phS = phSolver::Instance();

	// Get pointer to the single instance of PhGlobalArrayTransfer
	gat = PhGlobalArrayTransfer::Instance();

	// save the communicator
	iNewComm_C_ = MPI_COMM_WORLD;
	newcom.iNewComm = MPI_Comm_c2f(iNewComm_C_); // modifies newcom in fortran common block

	MPI_Comm_size(iNewComm_C_, &numProcs_);
	MPI_Comm_rank(iNewComm_C_, &rank_);

	// read configuration file
	input_fform();

	// update numProcsTotal and rank
//	MPI_Comm_size(newcomm, &numProcs);
//	MPI_Comm_rank(newcomm, &rank_);

	// Preprocess data and run the problem
	// Partition the problem to the correct number of processors

	if( rank_ == 0 )
	{
		//cout << "number of procs per simulation " << numprocs_perparticle << endl;
		Partition_Problem( numProcs_ );
	}

	MPI_Barrier(MPI_COMM_WORLD);

	sprintf(pathToProcsCaseDir,"%d-procs-case",numProcs_);
	chdir(pathToProcsCaseDir);

	input(&numProcs_, &rank_);
    proces();
    itrdrv_init(); // initialize solver

    // find nonzero entries of wall props
    for (int surfIdx = 0; surfIdx <= MAXSURF; surfIdx++)
    	if (nomodule.ValueListWallE[surfIdx] > 0)
    		WallEInd_.PushBack(surfIdx);

	// Compute the local state size
	// perhaps this can be less hard-coded
	// NSD is 3
	// 4 comes from the 3 vel components and 1 pressure component

	Nstate_local_ = (2 * 4) * conpar.nshguniq;

	if (nomodule.ideformwall > 0)
		Nstate_local_ += NSD * conpar.nshguniq;

	if (nomodule.imeasdist > 0)
		Nstate_local_ += conpar.nshguniq * 2; // distance vector and pre-multiplied distance vector

	if (rank_ == numProcs_ - 1) {

		Nstate_local_ += grcrbccom.numGRCRSrfs * 2; // number of lumped parameter states

		state_reduced_start_local_ = Nstate_local_; // index for first element of reduced state

		// Add the number of wall parameters to the local state size of the last process
		if (nreduced_has_wall_parameters_ && nomodule.ideformwall > 0)
			Nstate_local_ += nomodule.numWallRegions;

		if (nreduced_has_coupled_parameters_) {
			if (cp_rcr_estimate_compliance_)
				Nstate_local_ += grcrbccom.numGRCRSrfs;
			if (cp_rcr_estimate_resistance_)
				Nstate_local_ += grcrbccom.numGRCRSrfs;
		}
	}

	// Compute the reduced state size
	// (the number of variables to be estimated)
    if (nreduced_has_wall_parameters_)
    	Nreduced_ += nomodule.numWallRegions;

    if (nreduced_has_coupled_parameters_) {

    	if (cp_rcr_estimate_compliance_) {

    		int NCreduced = 0;

    		for (unsigned int kk = 0; kk < cp_rcr_include_compliance_.size(); kk++) {
    			if (cp_rcr_include_compliance_[kk])
    				NCreduced++;
    		}

    		//cout << "NCreduced " << NCreduced << endl;

    		//Nreduced_ += grcrbccom.numGRCRSrfs;
    		Nreduced_ += NCreduced;
    	}

    	if (cp_rcr_estimate_resistance_) {

    		int NRreduced = 0;

    		for (unsigned int kk = 0; kk < cp_rcr_include_resistance_.size(); kk++) {
    			if (cp_rcr_include_resistance_[kk])
    				NRreduced++;
            }

    		//cout << "NRreduced " << NRreduced << endl;

    		//Nreduced_ += grcrbccom.numGRCRSrfs;
    		Nreduced_ += NRreduced;
    	}
    }

	// Compute the global state size
	MPI_Allreduce(&Nstate_local_, &Nstate_, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	// Make sure the state array is distributed with the correct number
	// of elements per processor that is given by Nstate_local_
	duplicated_state_.Reallocate(Nstate_, Nstate_local_);

    param_out_.open ("estimated_reduced_state.dat");

    cout << "At rank " << rank_ << " we have " << Nstate_local_ << " state variables" << endl;

    if (rank_ == numProcs_ - 1) {
    	cout << "At rank " << rank_ << " we have " << Nreduced_ << " reduced variables" << endl;
    }
}


//! Initializes the first time step for the model.
void SimvascularVerdandiModel::InitializeFirstStep() {

}


//! Initializes the current time step for the model.
void SimvascularVerdandiModel::InitializeStep() {

}

void SimvascularVerdandiModel::Finalize() {

	itrdrv_finalize();

	param_out_.close();
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
void SimvascularVerdandiModel::ForwardFinalize() {

	itrdrv_iter_finalize(); // routines that allow moving to the next step

	// write down the parameter values in a file
	int state_start, state_end;

	duplicated_state_.GetProcessorRange(state_start, state_end);

	if (rank_ == numProcs_ -1) {

		int ncounter = 0;

		if (nreduced_has_wall_parameters_ && nomodule.ideformwall > 0) {
			for (int kk = 0; kk < nomodule.numWallRegions; kk++) {
				param_out_ << pow(2.0,duplicated_state_(state_start + state_reduced_start_local_ + ncounter) ) << " ";
				ncounter++;
			}
		}

		if (nreduced_has_coupled_parameters_ && grcrbccom.numGRCRSrfs > 0) {

			if (cp_rcr_estimate_compliance_)
				for (int kk = 0; kk < grcrbccom.numGRCRSrfs; kk++) {
					if (cp_rcr_include_compliance_[kk]) {
						param_out_ << pow(2.0,duplicated_state_(state_start + state_reduced_start_local_ + ncounter) ) << " ";
						ncounter++;
					}
				}

			if (cp_rcr_estimate_resistance_)
				for (int kk = 0; kk < grcrbccom.numGRCRSrfs; kk++) {
					if (cp_rcr_include_resistance_[kk]) {
						param_out_ << pow(2.0,duplicated_state_(state_start + state_reduced_start_local_ + ncounter) ) << " ";
						ncounter++;
					}
				}

		}

		param_out_ << endl;
	}
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

//	if (rank_ == numProcs_ -1)
//		cout << "begin particle step " << endl;

	duplicated_state_.Copy(x); // copies x into the duplicated state vector
	StateUpdated();            // updates actual model state with duplicated state vector


#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
	{
		stringstream s_temp;
		s_temp << this->rank_;
		string outname = "x_internal_before_apply-"+s_temp.str()+".dat";
		ofstream outfile(outname.c_str(), ofstream::app);
		outfile.precision(15);

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
		outfile << endl;
		local_X.Nullify();
	}
#endif

	itrdrv_iter_init();  // advances the model forward

	itrdrv_iter_step();  // note that ForwardFinalize is not called here

	x.Copy(GetState());        // copies the actual model state (via duplicated state) into x

//	if (rank_ == numProcs_ -1)
//		cout << "end particle step " << endl;

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
		cout << "getting state ";

	duplicated_state_.GetProcessorRange(state_start, state_end);

	icounter = state_start;

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

	if (nomodule.imeasdist > 0) {

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


	// get the states for the lumped parameter model

	if (rank_ == numProcs_ - 1) {

		for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++)
			duplicated_state_.SetBuffer(icounter++,(gat->global_lumped_parameter_Q)[surfIdx]);

		for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++)
			duplicated_state_.SetBuffer(icounter++,(gat->global_lumped_parameter_P)[surfIdx]);

	}

	// get the parameter value
	// note that this part of the state is only on the last processor

	if (rank_ == numProcs_ - 1) {

		if (nreduced_has_wall_parameters_ && nomodule.ideformwall > 0) {

			//val = nomodule.evw;
			//duplicated_state_.SetBuffer(icounter++,log2(val));

			for(int parIdx = 0; parIdx < nomodule.numWallRegions; parIdx++) {

				val = nomodule.ValueListWallE[WallEInd_(parIdx)];

				duplicated_state_.SetBuffer(icounter++,log2(val));

			}
			//cout << "[get] rank: " << rank_ << " val: " << log2(val) << endl;
		}

		if (nreduced_has_coupled_parameters_ && grcrbccom.numGRCRSrfs > 0) {

			if (cp_rcr_estimate_compliance_)
				for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++) {

					if (cp_rcr_include_compliance_[parIdx]) {
						val = gat->global_lumped_parameter_params[parIdx*3+1];

						duplicated_state_.SetBuffer(icounter++,log2(val));

						//				cout << "outlet " << parIdx << " ";
						//				cout << gat->global_lumped_parameter_params[parIdx*3+0] << " ";
						//			    cout << gat->global_lumped_parameter_params[parIdx*3+1] << " ";
						//				cout << gat->global_lumped_parameter_params[parIdx*3+2] << " " << endl;
					}

				}

			if (cp_rcr_estimate_resistance_)
				for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++) {

					if (cp_rcr_include_resistance_[parIdx]) {

						val = gat->global_lumped_parameter_params[parIdx*3+2];

						duplicated_state_.SetBuffer(icounter++,log2(val));

					}

				}

		}

	}

	duplicated_state_.Flush();

	if (rank_ == 0)
		cout << "[done]" << endl;

	return duplicated_state_;

}

//! Updates the internal state from the duplicated state
/*!

 */
void SimvascularVerdandiModel::StateUpdated() {

	int err;
	double val;
	int state_start, state_end, icounter;

	int actualIdx;

	if (rank_ == 0)
		cout << "setting state ";

	duplicated_state_.GetProcessorRange(state_start, state_end);

	icounter = state_start;

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

	}

	// set the P,Q states for the lumped parameter model

	if (rank_ == numProcs_ - 1) {

		// setting the Q state isn't strictly necessary
		// as Q is recomputed from the velocity field
		for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++) {

			(gat->global_lumped_parameter_Q)[surfIdx] = duplicated_state_(icounter++);

		}

		for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++) {

			(gat->global_lumped_parameter_P)[surfIdx] = duplicated_state_(icounter++);

		}
	}

	MPI_Bcast(gat->global_lumped_parameter_P,grcrbccom.numGRCRSrfs, MPI_DOUBLE, numProcs_-1, MPI_COMM_WORLD);
	MPI_Bcast(gat->global_lumped_parameter_Q,grcrbccom.numGRCRSrfs, MPI_DOUBLE, numProcs_-1, MPI_COMM_WORLD);
	//cout << "average distal P at rank " << rank_ << " : ";

//	for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++) {
//		cout << (gat->global_lumped_parameter_P)[surfIdx] << " ";
//	}
//	cout << endl;


	// set the values of the reduced state
	// note that this part of the state is only on the last processor
	// so we need to broadcast the parameter value to every processor

	if (nreduced_has_wall_parameters_ && nomodule.ideformwall > 0) {

		double* tempArray = new double[nomodule.numWallRegions];

		if (rank_ == numProcs_ - 1)
			for(int parIdx = 0; parIdx < nomodule.numWallRegions; parIdx++)
				tempArray[parIdx] = pow(2.0,duplicated_state_(icounter++));

		//MPI_Bcast(&val, 1, MPI_DOUBLE, numProcs_ - 1, MPI_COMM_WORLD);
		//nomodule.evw = val;

		MPI_Bcast(tempArray, nomodule.numWallRegions, MPI_DOUBLE, numProcs_ - 1, MPI_COMM_WORLD);

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

			MPI_Bcast(tempArray, grcrbccom.numGRCRSrfs, MPI_DOUBLE, numProcs_ - 1, MPI_COMM_WORLD);

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

			MPI_Bcast(tempArray, grcrbccom.numGRCRSrfs, MPI_DOUBLE, numProcs_ - 1, MPI_COMM_WORLD);

			for(int parIdx = 0; parIdx < grcrbccom.numGRCRSrfs; parIdx++)
				if (cp_rcr_include_resistance_[parIdx])
					gat->global_lumped_parameter_params[parIdx*3+2] = tempArray[parIdx];
		}

		delete [] tempArray;
	}

	// since we have set the state only on the master image nodes, we need to
	// copy state values from the master image to the interprocessor boundary nodes
	if (numProcs_ > 1)
		estim_helpers_setstate_comm();

    if (rank_ == 0)
		cout << "[done]" << endl;

}

////////////
// ERRORS //
////////////

/*! Returns a decomposition of the intial state error covariance matrix (\f$B\f$)
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

		if (nreduced_has_wall_parameters_)
			for (int i = 0; i < nomodule.numWallRegions; i++) {
				L.SetBuffer(start_ind_local + state_reduced_start_local_ + ncounter, ncounter, double(1));
				ncounter++;
			}

		if (nreduced_has_coupled_parameters_) {

			if(cp_rcr_estimate_compliance_)
				for (int i = 0; i < grcrbccom.numGRCRSrfs; i++) {
					if (cp_rcr_include_compliance_[i]) {
						L.SetBuffer(start_ind_local + state_reduced_start_local_ + ncounter, ncounter, double(1));
						ncounter++;
					}
				}

			if(cp_rcr_estimate_resistance_)
				for (int i = 0; i < grcrbccom.numGRCRSrfs; i++) {
					if (cp_rcr_include_resistance_[i]) {
						L.SetBuffer(start_ind_local + state_reduced_start_local_ + ncounter, ncounter, double(1));
						ncounter++;
					}
				}

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
		outfile.precision(15);

		U.WriteText(outfile);
	}
#endif


	// Read error variance value from a configuration file
	// Mlt(1./error_variance_value, U_);
//	L.Fill(0.);
//	for (int i = 0; i < Nreduced_; i++)
//		L(i + Nstate_ - Nreduced_, i) = 1.;

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
