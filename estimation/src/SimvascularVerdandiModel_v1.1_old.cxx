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
	//double saved_time = 0;
	//state saved_state;
	//if (!forward)
	//	saved_time = GetTime();

	//if (preserve_state)
	//	GetFullState(saved_state);

	this->SetState(x);
	phS->SolverForwardInit();
	phS->SolverForwardStep();
	this->GetState(x);   // for now, we will manually advance the time by calling the
	                     // the finalize step function of phsolver in the main
	                     // analysis loop

	//if (!forward)
	//	SetTime(saved_time);

	//if (preserve_state)
	//	SetFullState(saved_state);

	/*state x1(x), x2(x), x3(x);

	this->SetState(x1);
	phS->SolverForwardInit();
	phS->SolverForwardStep();
	this->GetState(x1);

	this->SetState(x2);
	phS->SolverForwardInit();
	phS->SolverForwardStep();
	this->GetState(x2);

	this->SetState(x3);
	phS->SolverForwardInit();
	phS->SolverForwardStep();
	this->GetState(x3);

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

	this->SetState(x);
	phS->SolverForwardInit();
	phS->SolverForwardStep();
	this->GetState(x);*/
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
	return (double)phS->getTime();
}


//! Sets the time of the model to a given time.
/*!
      \param[in] time a given time.
 */
void SimvascularVerdandiModel::SetTime(double time)
{
	phS->setTime((int)time);
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


//! Provides the state vector.
/*!
      \param[out] state the reduced state vector.
      this function takes values from the variables yold,acold,uold
 */
void SimvascularVerdandiModel::GetState(state& state)
{
	// transfer the solution field
	//SCField solField;
	int err,icounter;
	double val;

	cout << "getting state" << endl;

	state.Reallocate(Nstate_);

	icounter = 0;

	// get the "solution field" we should probably do this more efficiently than copying
	// in the flowsolver this is yold
	for(int unitIdx=0; unitIdx < phS->GetRequiredField(0)->GetNumUnits(); unitIdx++) {
		for(int varIdx=0; varIdx < phS->GetRequiredField(0)->GetNumVars()-1; varIdx++) { // ignore the 5th dof and beyond
			phS->GetValueDbl(*phS->GetRequiredField(0), unitIdx, varIdx, val);

			state(icounter) = val;

			icounter++;
		}
	}

	// get the "acceleration field"
	// in the flowsolver this is acold
	for(int unitIdx=0; unitIdx < phS->GetRequiredField(1)->GetNumUnits(); unitIdx++) {
		for(int varIdx=0; varIdx < phS->GetRequiredField(1)->GetNumVars()-1; varIdx++) { // ignore the 5th dof and beyond
			phS->GetValueDbl(*phS->GetRequiredField(1), unitIdx, varIdx, val);

			state(icounter) = val;

			icounter++;
		}
	}

	// get the "displacement field"
	// in the flowsolver this is uold
	for(int unitIdx=0; unitIdx < phS->GetRequiredField(2)->GetNumUnits(); unitIdx++) {
		for(int varIdx=0; varIdx < phS->GetRequiredField(2)->GetNumVars(); varIdx++) { // 3 dofs
			phS->GetValueDbl(*phS->GetRequiredField(2), unitIdx, varIdx, val);

			state(icounter) = val;

			icounter++;
		}

	}

	// get RCR history ?

	// get lagrange multipliers ?

	// get the elastic modulus and take the log of it -- we don't want to restrict the parameter to only take
	// positive values
	// in the flowsolver this is evw
	phS->GetValueDbl(*phS->GetRequiredField(3), 0, 0, val);
	state(icounter) = log(val); // temporary

	cout << "done" << endl;

}

//! Sets the state vector.
/*! Before setting the reduced state vector, special requirements can be
      enforced; e.g. positivity requirement or inferior and superior limits.
      \param[in] state the reduced state vector.
 */
void SimvascularVerdandiModel::SetState(state& state)
{
	int icounter = 0;
	double val;

	cout << "setting state" << endl;

	// set the "solution field" we should probably do this more efficiently than copying
	for(int unitIdx=0; unitIdx < phS->GetRequiredField(0)->GetNumUnits(); unitIdx++) {
		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond

			val = state(icounter);

			phS->SetValueDbl(*phS->GetRequiredField(0), unitIdx, varIdx, val);

			icounter++;
		}
	}

	// set the "acceleration field"
	for(int unitIdx=0; unitIdx < phS->GetRequiredField(1)->GetNumUnits(); unitIdx++) {
		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond

			val = state(icounter);

			phS->SetValueDbl(*phS->GetRequiredField(1), unitIdx, varIdx, val);

			icounter++;
		}
	}

	// set the "displacement field"
	for(int unitIdx=0; unitIdx < phS->GetRequiredField(2)->GetNumUnits(); unitIdx++) {
		for(int varIdx=0; varIdx < phS->GetRequiredField(2)->GetNumVars(); varIdx++) { // 3 dofs

			val = state(icounter);

			phS->SetValueDbl(*phS->GetRequiredField(2), unitIdx, varIdx, val);

			icounter++;
		}

	}

	// set RCR history ?

	// set lagrange multipliers ?

	// set the elastic modulus by taking the exp of the parameter

	val = exp(state(icounter));
	phS->SetValueDbl(*phS->GetRequiredField(3), 0, 0, val);

	cout << "done" << endl;
}

//! Provides the full state vector.
/*!
      \param[out] state the full state vector.
 */
void SimvascularVerdandiModel::GetFullState(state& state)
{
	GetState(state);
}


//! Sets the full state vector.
/*!
      \param[in] state the full state vector.
 */
void SimvascularVerdandiModel::SetFullState(state& state)
{
	SetState(state);
}

////////////
// ERRORS //
////////////

/*! Returns a decomposition of the state error covariance matrix (\f$B\f$)
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
	U.SetIdentity();
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
