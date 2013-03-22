#ifndef VERDANDI_FILE_MODEL_SimvascularVerdandiModel_CXX

#include "SimvascularVerdandiModel.hxx"

namespace Verdandi {


////////////////////////////////
// CONSTRUCTOR AND DESTRUCTOR //
////////////////////////////////


//! Constructor.
SimvascularVerdandiModel::SimvascularVerdandiModel() {
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
	configuration.Set("error_statistics.state_error_variance",state_error_variance_value_);
	//cout << "assimilation file: " << configuration.GetFilePath() << endl;

    Initialize();

}

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
	gat = GlobalArrayTransfer::Instance();

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

	// Initialize Nstate_, Nparameter_, state_....
	Nparameter_ = 1;

	// Compute the local state size
	// perhaps this can be less hard-coded
	// NSD is 3
	// 4 comes from the 3 vel components and 1 pressure component
	Nstate_local_ = (NSD + 2 * 4) * conpar.nshguniq;

	if (rank_ == numProcs_ - 1) {
		// Add the number of lumped parameter surfaces (twice) to the local state size of the last process
		Nstate_local_ += grcrbccom.numGRCRSrfs * 2;

		// Add the number of parameters to the local state size of the last process
		Nstate_local_ += Nparameter_;
	}

	// Compute the global state size
	MPI_Allreduce(&Nstate_local_, &Nstate_, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	// Make sure the state array is distributed with the correct number
	// of elements per processor that is given by Nstate_local_
	//
	duplicated_state_.Reallocate(Nstate_, Nstate_local_);

    param_out_.open ("estimated_parameters_mean.dat");

    if (rank_ == 0)
    	cout << "Simvascular Model initiated" << endl;
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

	//
	// write down the parameter values in a file
	//

	if (rank_ == numProcs_ -1) {
		if (nomodule.ideformwall > 0) {
			for (int kk = 0; kk < Nparameter_; kk++) {
				param_out_ << pow(2.0,duplicated_state_(Nstate_-Nparameter_+kk)) << " ";
			}
			param_out_ << endl;
		}
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

	duplicated_state_.Copy(x); // copies x into the duplicated state vector
	StateUpdated();            // updates actual model state with duplicated state vector

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

	MPI_Barrier(MPI_COMM_WORLD);

	duplicated_state_.GetProcessorRange(state_start, state_end);

	icounter = state_start;

	//
	// Get the "solution field"
	// The state elements must be distributed such that there are no duplicate entries across processors
	// The list of unique nodes (master image) is subset of the nodes that are on the local processor
	//

	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

		actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond

			duplicated_state_.SetBuffer(icounter++,(gat->global_yold_ptr)[varIdx * conpar.nshg + actualIdx-1]);

		}

	}

	//
	// get the "acceleration field"
	// in the flowsolver this is acold
	//

	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

		actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond

			duplicated_state_.SetBuffer(icounter++,(gat->global_acold_ptr)[varIdx * conpar.nshg + actualIdx-1]);

		}

	}

	//
	// get the "displacement field"
	// in the flowsolver this is uold
	//

	if (nomodule.ideformwall > 0) {

		for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

			actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

			for(int varIdx=0; varIdx < 3; varIdx++) { // 3 dofs

				duplicated_state_.SetBuffer(icounter++,(gat->global_uold_ptr)[varIdx * conpar.nshg + actualIdx-1]);

			}

		}

	}


	//
	// get the P,Q states for the lumped parameter model
	//

	if (rank_ == numProcs_ - 1) {

		for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++) {

			duplicated_state_.SetBuffer(icounter++,(gat->global_lumped_parameter_P)[surfIdx]);

		}

		for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++) {

			duplicated_state_.SetBuffer(icounter++,(gat->global_lumped_parameter_Q)[surfIdx]);

		}

	}

	//
	// get the parameter value
	// note that this part of the state is only on the last processor
	//

	if (rank_ == numProcs_ - 1) {

		if (nomodule.ideformwall > 0) {

			val = nomodule.evw;

			duplicated_state_.SetBuffer(icounter++,log2(val));

			//cout << "[get] rank: " << rank_ << " val: " << log2(val) << endl;
		}

	}

	duplicated_state_.Flush();

	if (rank_ == 0)
		cout << "[done]" << endl;

//	duplicated_state_.Print();

	MPI_Barrier(MPI_COMM_WORLD);

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

	MPI_Barrier(MPI_COMM_WORLD);

	duplicated_state_.GetProcessorRange(state_start, state_end);

	icounter = state_start;

	//
	// set the "solution field"
	// The list of unique nodes (master image) is subset of the nodes that are on the local processor
	//

	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

		actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond

			(gat->global_yold_ptr)[varIdx * conpar.nshg + actualIdx-1] = duplicated_state_(icounter++);

		}

	}

	//
	// set the "acceleration field"
	// in the flowsolver this is acold
	//

	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

		actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond

			(gat->global_acold_ptr)[varIdx * conpar.nshg + actualIdx-1] = duplicated_state_(icounter++);

		}

	}

	//
	// set the "displacement field"
	// in the flowsolver this is uold
	//

	if (nomodule.ideformwall > 0) {

		for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

			actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

			for(int varIdx=0; varIdx < 3; varIdx++) { // 3 dofs

				(gat->global_uold_ptr)[varIdx * conpar.nshg + actualIdx-1] = duplicated_state_(icounter++);

			}

		}

	}

	//
	// set the P,Q states for the lumped parameter model
	//

	if (rank_ == numProcs_ - 1) {

		for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++) {

			(gat->global_lumped_parameter_P)[surfIdx] = duplicated_state_(icounter++);

		}

		for(int surfIdx=0; surfIdx < grcrbccom.numGRCRSrfs; surfIdx++) {

			(gat->global_lumped_parameter_Q)[surfIdx] = duplicated_state_(icounter++);

		}

	}

	//
	// set the parameter value
	// note that this part of the state is only on the last processor
	// so we need to broadcast the parameter value to every processor
	//

	if (rank_ == numProcs_ - 1) {

		if (nomodule.ideformwall > 0) {
			val = pow(2.0,duplicated_state_(icounter));
		}

		//cout << "[set] rank: " << rank_ << " val: " << duplicated_state_(icounter) << endl;
	}

	if (nomodule.ideformwall > 0) {
		MPI_Bcast(&val, 1, MPI_DOUBLE, numProcs_ - 1, MPI_COMM_WORLD);

		nomodule.evw = val;

	}

	//
	// since we have set the state only on the master image nodes, we need to
	// copy state values from the master image to the interprocessor boundary nodes
	//

	if (numProcs_ > 1) {

		estim_helpers_setstate_comm();

	}

	//
	// lumped parameter model states
	//

    if (rank_ == 0)
		cout << "[done]" << endl;

    MPI_Barrier(MPI_COMM_WORLD);

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

	//
	// L is distributed across processor (since one of its dimensions is Nstate)
	// but U is on every processor (dimensions of Nparameter)
	//

	L.Reallocate(Nstate_, Nparameter_, Nstate_local_);

	if (rank_ == numProcs_ - 1)  // it's zero except in the last processor
		for (int i = 0; i < Nparameter_; i++)
			L.SetBuffer(i + Nstate_ - Nparameter_, i, double(1));
	L.Flush();

	U.Reallocate(Nparameter_, Nparameter_);
	U.Fill(double(0));

	//state_error_variance_value_ = double(1);

	for (int i = 0; i < Nparameter_; i++)
		U(i, i) = double(double(1) / state_error_variance_value_);

	// Read error variance value from a configuration file
	// Mlt(1./error_variance_value, U_);
//	L.Fill(0.);
//	for (int i = 0; i < Nparameter_; i++)
//		L(i + Nstate_ - Nparameter_, i) = 1.;

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
