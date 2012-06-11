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

	MPI_Finalize();
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
	int rank;
	int numProcsTotal,numProcs;
	int ierr = 0;
	char inpfilename[100];

	MPI_Comm newcomm;

	int color,key;

	int numparticles = 1;
	int numprocs_perparticle;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcsTotal);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// create MPI communicator for one simulation
	// by splitting MPI_COMM_WORLD
	numprocs_perparticle = numProcsTotal / numparticles;
	color = rank / numprocs_perparticle;
	key = 0;
	MPI_Comm_split(MPI_COMM_WORLD,color,key,&newcomm);

	// Initialize phSolver
	phS = phSolver::Instance();

	// save the communicator
	phS->setCommunicator(newcomm);

	// read configuration file
	phS->readConfiguration();

	// update numProcsTotal and rank
	MPI_Comm_size(newcomm, &numProcs);
	MPI_Comm_rank(newcomm, &rank);

	// Preprocess data and run the problem
	// Partition the problem to the correct number of processors

	if( rank == 0 )
	{
		cout << "number of procs per simulation " << numprocs_perparticle << endl;
		Partition_Problem( numProcsTotal, color, numprocs_perparticle );
	}

	MPI_Barrier(MPI_COMM_WORLD);

	sprintf(inpfilename,"%d-procs-%d-simid",numprocs_perparticle,color);

	cout << "changing directory to " << inpfilename << endl;

	phS->readMeshAndSolution_fromFiles(inpfilename);

	solveReturn = phS->SolverInit(); // initialize solver

	if ( 0 ==solveReturn ) {
		if (rank == 0)
			cout << "flowsolver initiated" << endl;
	}
	else {
		if (rank == 0)
			fprintf(stderr, "flowsolver failed to initiate\n");
		//ierr = 1;
		MPI_Comm_free(&newcomm);
		MPI_Finalize();
		exit(1);
	}

	// Initialize Nstate_, Nparameter_, state_....
	Nparameter_ = 1;
	Nstate_ = phS->getSize()+Nparameter_;
	Nobservation_ = 6;

	duplicated_state_.Reallocate(Nstate_);

	//x_.Reallocate(Nstate_);

	//Matrix<double> H(Nstate_, Nobservation_);
	//H.SetIdentity();

	//H.Write("H.bin");

	//Matrix<double, General, ArrayRowSparse> H(Nstate_, Nobservation_);
	//H.AddValue(i, j, value);
	//H.Write
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
	phS->SolverForwardFinalize();
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

	if (preserve_state)
		saved_state.SetData(duplicated_state_);

	duplicated_state_.Nullify();
	duplicated_state_.SetData(x);
	StateUpdated();

	phS->SolverForwardInit();
	phS->SolverForwardStep();

	GetStateCopy (duplicated_state_);
	duplicated_state_.Nullify();

	// For now the "forward" functionality is partially implemented in "ForwardFinalize"
	// which should be called in the main loop of the ROUKF driver

	//if (!forward)
	//	SetTime(saved_time);

	if (preserve_state) {
		duplicated_state_.SetData(saved_state);
		StateUpdated();
		saved_state.Nullify();
	}

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


//! Returns the size of the full state vector.
/*!
      \return The size of the full state vector.
 */
int SimvascularVerdandiModel::GetNfull_state() const
{
	return GetNstate();
}


//! copies internal state to "state"
/*!
      \param[in] state: a reference to state vector
 */
void SimvascularVerdandiModel::GetStateCopy(state& state)
{
	int err,icounter;
	double val;

	cout << "getting state ";

	state.Reallocate(Nstate_);

	icounter = 0;

	// get the "solution field" we should probably do this more efficiently than copying
	// in the flowsolver this is yold
	for(int unitIdx=0; unitIdx < phS->GetProducedField(0)->GetNumUnits(); unitIdx++) {
		for(int varIdx=0; varIdx < phS->GetProducedField(0)->GetNumVars()-1; varIdx++) { // ignore the 5th dof and beyond
			phS->GetValueDbl(*phS->GetProducedField(0), unitIdx, varIdx, val);

			state(icounter) = val;

			icounter++;
		}
	}

	// get the "acceleration field"
	// in the flowsolver this is acold
	for(int unitIdx=0; unitIdx < phS->GetProducedField(1)->GetNumUnits(); unitIdx++) {
		for(int varIdx=0; varIdx < phS->GetProducedField(1)->GetNumVars()-1; varIdx++) { // ignore the 5th dof and beyond
			phS->GetValueDbl(*phS->GetProducedField(1), unitIdx, varIdx, val);

			state(icounter) = val;

			icounter++;
		}
	}

	// get the "displacement field"
	// in the flowsolver this is uold
	for(int unitIdx=0; unitIdx < phS->GetProducedField(2)->GetNumUnits(); unitIdx++) {
		for(int varIdx=0; varIdx < phS->GetProducedField(2)->GetNumVars(); varIdx++) { // 3 dofs
			phS->GetValueDbl(*phS->GetProducedField(2), unitIdx, varIdx, val);

			state(icounter) = val;

			icounter++;
		}

	}

	// get RCR history ?

	// get lagrange multipliers ?

	// get the elastic modulus and take the log of it -- we don't want to restrict the parameter to only take
	// positive values
	// in the flowsolver this is evw
	phS->GetValueDbl(*phS->GetProducedField(3), 0, 0, val);
	state(icounter) = log2(val); // temporary

	cout << "[done]" << endl;

}

//! Sets the internal state from "state".
/*!
      \param[in] state: reference to state vector.
 */
void SimvascularVerdandiModel::SetStateCopy(state& state)
{
	int icounter = 0;
	double val;

	cout << "setting state ";

	// set the "solution field" we should probably do this more efficiently than copying
	for(int unitIdx=0; unitIdx < phS->GetProducedField(0)->GetNumUnits(); unitIdx++) {
		for(int varIdx=0; varIdx < phS->GetProducedField(0)->GetNumVars()-1; varIdx++) { // ignore the 5th dof and beyond

			val = state(icounter);

			phS->SetValueDbl(*phS->GetProducedField(0), unitIdx, varIdx, val);

			icounter++;
		}
	}

	// set the "acceleration field"
	for(int unitIdx=0; unitIdx < phS->GetProducedField(1)->GetNumUnits(); unitIdx++) {
		for(int varIdx=0; varIdx < phS->GetProducedField(1)->GetNumVars()-1; varIdx++) { // ignore the 5th dof and beyond

			val = state(icounter);

			phS->SetValueDbl(*phS->GetProducedField(1), unitIdx, varIdx, val);

			icounter++;
		}
	}

	// set the "displacement field"
	for(int unitIdx=0; unitIdx < phS->GetProducedField(2)->GetNumUnits(); unitIdx++) {
		for(int varIdx=0; varIdx < phS->GetProducedField(2)->GetNumVars(); varIdx++) { // 3 dofs

			val = state(icounter);

			phS->SetValueDbl(*phS->GetProducedField(2), unitIdx, varIdx, val);

			icounter++;
		}

	}

	// set RCR history ?

	// set lagrange multipliers ?

	// set the elastic modulus by taking the exp of the parameter

	val = pow(2.0,state(icounter));
	phS->SetValueDbl(*phS->GetProducedField(3), 0, 0, val);

	cout << "[done]" << endl;
}

//! Provides the full state vector.
/*!
      \param[in] state: reference to state vector
 */
void SimvascularVerdandiModel::GetFullStateCopy(state& state)
{

}


//! Sets the full state vector.
/*!
      \param[in] state: reference to state vector
 */
void SimvascularVerdandiModel::SetFullStateCopy(state& state)
{
	//SetStateCopy(state);
}

//! Returns a reference to the duplicated state.
/*!
      \return refernce a to duplicated state
 */
SimvascularVerdandiModel::state& SimvascularVerdandiModel::GetState()
{
    GetStateCopy(duplicated_state_);

	return duplicated_state_;
}

//! Updates the internal state from the duplicated state
/*!

 */
void SimvascularVerdandiModel::StateUpdated()
{
    SetStateCopy(duplicated_state_);
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
void SimvascularVerdandiModel::GetStateErrorVarianceSqrt(state_error_variance& L,
		state_error_variance& U)
{
	L.Reallocate(Nstate_, Nparameter_), U.Reallocate(Nparameter_, Nparameter_);
	U.SetIdentity(); // the inverse of the parameter covariance matrix
	// Read error variance value from a configuration file
	// Mlt(1./error_variance_value, U_);
	L.Fill(0.);
	for (int i = 0; i < Nparameter_; i++)
		L(i + Nstate_ - Nparameter_, i) = 1.;

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
