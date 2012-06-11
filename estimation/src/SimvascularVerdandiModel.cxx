#ifndef VERDANDI_FILE_MODEL_SimvascularVerdandiModel_CXX

#include "SimvascularVerdandiModel.hxx"

namespace Verdandi
{


////////////////////////////////
// CONSTRUCTOR AND DESTRUCTOR //
////////////////////////////////


//! Constructor.
SimvascularVerdandiModel::SimvascularVerdandiModel()
{
}


//! Destructor.
SimvascularVerdandiModel::~SimvascularVerdandiModel()
{
	phS->SolverFinalize();

	param_out_.close();

	//MPI_Finalize();
}


////////////////
// INITIALIZE //
////////////////


//! Initializes the model.
/*!
      \param[in] configuration_file configuration file.
 */
void SimvascularVerdandiModel::Initialize(string configuration_file)
{
	throw ErrorUndefined("SimvascularVerdandiModel"
			"::Initialize(string configuration_file)");
}

void SimvascularVerdandiModel::Initialize(int argc, char * argv[])
{
	char inpfilename[100];

//	MPI_Comm newcomm;
//	int color,key;
//	int numparticles = 1;
//	int numprocs_perparticle;

	//MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs_);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

	// create MPI communicator for one simulation
	// by splitting MPI_COMM_WORLD
	// this functionality isn't used for now as
	// the "numparticles" always is 1
//	numprocs_perparticle = numProcsTotal / numparticles;
//	color = rank_ / numprocs_perparticle;
//	key = 0;
//	MPI_Comm_split(MPI_COMM_WORLD,color,key,&newcomm);

	// Initialize phSolver
	phS = phSolver::Instance();

	// save the communicator
	phS->setCommunicator(MPI_COMM_WORLD);

	// read configuration file
	phS->readConfiguration();

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

	sprintf(inpfilename,"%d-procs-case",numProcs_);

	cout << "changing directory to " << inpfilename << endl;

	phS->readMeshAndSolution_fromFiles(inpfilename);

	solveReturn = phS->SolverInit(); // initialize solver

	if ( 0 ==solveReturn ) {
		if (rank_ == 0)
			cout << "flowsolver initiated" << endl;
	}
	else {
		if (rank_ == 0)
			fprintf(stderr, "flowsolver failed to initiate\n");
		MPI_Finalize();
		exit(1);
	}

	// Initialize Nstate_, Nparameter_, state_....
	Nparameter_ = 1;

	// Compute the local state size
	Nstate_local_ = phS->getSize();

	// Add the number of parameters to the local state size of the last processor
	if (rank_ == numProcs_ - 1)
		Nstate_local_ += Nparameter_;

	// Compute the global state size
	MPI_Allreduce(&Nstate_local_, &Nstate_, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	// Make sure the state array is distributed with the correct number
	// of elements per processor that is given by Nstate_local_
	//
	duplicated_state_.Reallocate(Nstate_, Nstate_local_);

	node_field_ = phS->GetRequiredField("local index of unique nodes");
    soln_field_ = phS->GetRequiredField("solution");
    acc_field_ = phS->GetRequiredField("time derivative of solution");
    if (nomodule.ideformwall > 0)
    	disp_field_ = phS->GetRequiredField("displacement");

	//x_.Reallocate(Nstate_);

	//Matrix<double> H(Nstate_, Nobservation_);
	//H.SetIdentity();

	//H.Write("H.bin");

	//Matrix<double, General, ArrayRowSparse> H(Nstate_, Nobservation_);
	//H.AddValue(i, j, value);
	//H.Write

    param_out_.open ("estimated_parameters_mean.dat");

    if (rank_ == 0)
    	cout << "Simvascular Model initiated" << endl;


//    int tempcount = 0,tempval = 0;
//    for (int kk = 0; kk < phS->GetNumBlocks(); kk++) {
//    	tempcount += phS->GetBlockSize(kk);
//    	cout << tempcount << endl;
//    }
//    cout << "Number of blocks: " << phS->GetNumBlocks() << " Number of elements total: " << tempcount << endl;
}


//! Initializes the first time step for the model.
void SimvascularVerdandiModel::InitializeFirstStep()
{

}


//! Initializes the current time step for the model.
void SimvascularVerdandiModel::InitializeStep()
{

}


////////////////
// PROCESSING //
////////////////


//! Advances one step forward in time.
/*! \f[x^f_{h+1} = \mathcal{M}_h(x^a_h, p_h)\,.\f] */
void SimvascularVerdandiModel::Forward()
{
	phS->SolverForwardInit();
	phS->SolverForwardStep();
	phS->SolverForwardFinalize();
}

//! Finalizes the current time step
//! meant to be called after the mean state with innovation is set
void SimvascularVerdandiModel::ForwardFinalize()
{
	phS->SolverForwardFinalize(); // routines that allow moving to the next step

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
bool SimvascularVerdandiModel::HasFinished() const
{
	return phS->hasFinished();
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
		bool forward, bool preserve_state)
{
	double saved_time = 0;
	state saved_state;
	//if (!forward)
	//	saved_time = GetTime();

//	if (preserve_state)
//		GetStateCopy(saved_state);

	duplicated_state_.Copy(x); // copies x into the duplicated state vector
	StateUpdated();            // updates actual model state with duplicated state vector

	phS->SolverForwardInit();  // advances the model forward
	phS->SolverForwardStep();  // note that ForwardFinalize is not called here

	x.Copy(GetState());        // copies the actual model state (via duplicated state) into x

	// For now the "forward" functionality is partially implemented in "ForwardFinalize"
	// which should be called in the main loop of the ROUKF driver

	//if (!forward)
	//	SetTime(saved_time);

//	if (preserve_state)
//	{
//		state_.Copy(saved_state);
//		StateUpdated();
//	}

	/*state x1(x), x2(x), x3(x);

	this->SetStateCopy(x1);
	phS->SolverForwardInit();
	phS->SolverForwardStep();
	this->GetStateCopy(x1);

	this->SetStateCopy(x2);
	phS->SolverForwardInit();
	phS->SolverForwardStep();
	this->GetStateCopy(x2);

	this->SetStateCopy(x3);
	phS->SolverForwardInit();
	phS->SolverForwardStep();
	this->GetStateCopy(x3);

	for (int i = 0; i < x1.GetM(); i++)
		if (x1(i) != x2(i))
			throw ErrorProcessing("CheckingModel<Model>::ApplyOperator"
					"(state& x,bool forward,"
					" bool preserve_state)",
					"x1 = x2 but ApplyOperator(x1) != "
							"ApplyOperator(x2)\n"
							"x1 = x2 = " + to_str(x1) + "\n"
							"ApplyOperator(x1)(" + to_str(i) + ") = "
							+ to_str(x1(i)) + "\n ApplyOperator(x2)("
							+ to_str(i) + ") = " + to_str(x2(i)) + ".");
	for (int i = 0; i < x1.GetM(); i++)
		if (x1(i) != x3(i))
			throw ErrorProcessing("CheckingModel<Model>::ApplyOperator"
					"(state& x,bool forward,"
					" bool preserve_state)",
					"x1 = x3 but ApplyOperator(x1) != "
							"ApplyOperator(x3)\n"
							"x1 = x3 = " + to_str(x1) + "\n"
							"ApplyOperator(x1)(" + to_str(i) + ") = "
							+ to_str(x1(i)) + "\n ApplyOperator(x2)("
							+ to_str(i) + ") = " + to_str(x3(i)) + ".");

	this->SetStateCopy(x);
	phS->SolverForwardInit();
	phS->SolverForwardStep();
	this->GetStateCopy(x);*/
}


////////////////////
// ACCESS METHODS //
////////////////////


//! Returns the current time.
/*!
      \return The current time.
 */
double SimvascularVerdandiModel::GetTime() const
{
	// the time is adjusted by 1 due to the time not being
	// incremented until the very end of the time step
	return (double)(phS->getTime()+1);
}


//! Sets the time of the model to a given time.
/*!
      \param[in] time a given time.
 */
void SimvascularVerdandiModel::SetTime(double time)
{
	throw ErrorUndefined("SimvascularVerdandiModel"
			"::SetTime(double time)");
}


//! Returns the state vector size.
/*!
      \return The state vector size.
 */
int SimvascularVerdandiModel::GetNstate() const
{
	return Nstate_;
}


//! Returns the size of the local (on-processor) state vector.
/*!
      \return The size of the local state vector.
 */
int SimvascularVerdandiModel::GetLocalNstate() const
{
	return Nstate_local_; // this is used in reallocate routine for petsc matrices in ROUKF
}


//! copies internal state to "duplicated state" and returns a reference to it
/*!
      \return: a reference to "duplicated state" vector
 */
SimvascularVerdandiModel::state& SimvascularVerdandiModel::GetState()
{
	int err;
	double val;
	int state_start, state_end, icounter;

	int actualIdx;

	if (rank_ == 0)
		cout << "getting state ";

	duplicated_state_.GetProcessorRange(state_start, state_end);

	icounter = state_start;

	//
	// Get the "solution field"
	// The state elements must be distributed such that there are no duplicate entries across processors
	// The list of unique nodes (master image) is subset of the nodes that are on the local processor
	//

	for(int unitIdx=0; unitIdx < node_field_->GetNumUnits(); unitIdx++) {

		phS->GetValue(*node_field_, unitIdx, 0, actualIdx);

		for(int varIdx=0; varIdx < soln_field_->GetNumVars()-1; varIdx++) { // ignore the 5th dof and beyond

			phS->GetValue(*soln_field_, actualIdx-1, varIdx, val);

			duplicated_state_.SetBuffer(icounter++,val);

		}

	}

	//
	// get the "acceleration field"
	// in the flowsolver this is acold
	//

	for(int unitIdx=0; unitIdx < node_field_->GetNumUnits(); unitIdx++) {

		phS->GetValue(*node_field_, unitIdx, 0, actualIdx);

		for(int varIdx=0; varIdx < acc_field_->GetNumVars()-1; varIdx++) { // ignore the 5th dof and beyond

			phS->GetValue(*acc_field_, actualIdx-1, varIdx, val);

			duplicated_state_.SetBuffer(icounter++,val);

		}

	}

	//
	// get the "displacement field"
	// in the flowsolver this is uold
	//

	if (nomodule.ideformwall > 0) {

		for(int unitIdx=0; unitIdx < node_field_->GetNumUnits(); unitIdx++) {

			phS->GetValue(*node_field_, unitIdx, 0, actualIdx);

			for(int varIdx=0; varIdx < disp_field_->GetNumVars(); varIdx++) { // 3 dofs

				phS->GetValue(*disp_field_, actualIdx-1, varIdx, val);

				duplicated_state_.SetBuffer(icounter++,val);

			}

		}

	}

	//
	// get the parameter value
	// note that this part of the state is only on the last processor
	//

	if (rank_ == numProcs_ - 1) {

		if (nomodule.ideformwall > 0) {
			phS->GetValue(*phS->GetRequiredField("elastic modulus scalar"), 0, 0, val);

			duplicated_state_.SetBuffer(icounter++,log2(val));
		}

	}

	duplicated_state_.Flush();

	if (rank_ == 0)
		cout << "[done]" << endl;

//	duplicated_state_.Print();

	return duplicated_state_;

}

//! Updates the internal state from the duplicated state
/*!

 */
void SimvascularVerdandiModel::StateUpdated()
{
	int err;
	double val;
	int state_start, state_end, icounter;

	int actualIdx;

	if (rank_ == 0)
		cout << "setting state ";

	duplicated_state_.GetProcessorRange(state_start, state_end);

//	duplicated_state_.Print();

	icounter = state_start;

	//
	// set the "solution field"
	// The list of unique nodes (master image) is subset of the nodes that are on the local processor
	//

	for(int unitIdx=0; unitIdx < node_field_->GetNumUnits(); unitIdx++) {

		phS->GetValue(*node_field_, unitIdx, 0, actualIdx);

		for(int varIdx=0; varIdx < soln_field_->GetNumVars()-1; varIdx++) { // ignore the 5th dof and beyond

			val = duplicated_state_(icounter++);

			phS->SetValue(*soln_field_, actualIdx-1, varIdx, val);

		}

	}

	//
	// set the "acceleration field"
	// in the flowsolver this is acold
	//

	for(int unitIdx=0; unitIdx < node_field_->GetNumUnits(); unitIdx++) {

		phS->GetValue(*node_field_, unitIdx, 0, actualIdx);

		for(int varIdx=0; varIdx < acc_field_->GetNumVars()-1; varIdx++) { // ignore the 5th dof and beyond

			val = duplicated_state_(icounter++);

			phS->SetValue(*acc_field_, actualIdx-1, varIdx, val);

		}

	}

	//
	// set the "displacement field"
	// in the flowsolver this is uold
	//

	for(int unitIdx=0; unitIdx < node_field_->GetNumUnits(); unitIdx++) {

		phS->GetValue(*node_field_, unitIdx, 0, actualIdx);

		for(int varIdx=0; varIdx < disp_field_->GetNumVars(); varIdx++) { // 3 dofs

			val = duplicated_state_(icounter++);

			phS->SetValue(*disp_field_, actualIdx-1, varIdx, val);

		}

	}

	//
	// set the parameter value
	// note that this part of the state is only on the last processor
	// so we need to broadcast the parameter value to every processor
	//

	if (rank_ == numProcs_ - 1) {

		if (nomodule.ideformwall > 0) {
			val = pow(2.0,duplicated_state_(icounter++));
		}

	}

	if (nomodule.ideformwall > 0) {
		MPI_Bcast(&val, 1, MPI_DOUBLE, numProcs_ - 1, MPI_COMM_WORLD);

		phS->SetValue(*phS->GetRequiredField("elastic modulus scalar"), 0, 0, val);
	}

	//
	// since we have set the state only on the master image nodes, we need to
	// copy state values from the master image to the interprocessor boundary nodes
	//

	if (numProcs_ > 1) {

		setstate_comm();

	}

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
void SimvascularVerdandiModel::GetStateErrorVarianceSqrt(L_matrix& L, U_matrix& U)
{

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

	state_error_variance_value_ = double(1);

	for (int i = 0; i < Nparameter_; i++)
		U(i, i) = double(double(1) / state_error_variance_value_);

	// Read error variance value from a configuration file
	// Mlt(1./error_variance_value, U_);
//	L.Fill(0.);
//	for (int i = 0; i < Nparameter_; i++)
//		L(i + Nstate_ - Nparameter_, i) = 1.;

}


//! Returns the name of the class.
/*!
      \return The name of the class.
 */
string SimvascularVerdandiModel::GetName() const
{
	return "SimvascularVerdandiModel";
}


//! Receives and handles a message.
/*
      \param[in] message the received message.
 */
void SimvascularVerdandiModel::Message(string message)
{
	// Put here any processing you need.
}


}

#define VERDANDI_FILE_MODEL_SimvascularVerdandiModel_CXX
#endif
