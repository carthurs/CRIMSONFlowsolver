#ifndef VERDANDI_FILE_MODEL_SimvascularVerdandiModel_HXX

#include "common_c.h"yy

#include "partition.h"
#include "input_fform.h"
#include "input.h"
#include "proces.h"
#include "itrdrv.h"
#include "estimation_helpers.h"
#include "GlobalArrayTransfer.h"

#include "mpi.h"

#ifdef intel
#include <direct.h>
#else
#include <unistd.h>
#endif

#include <iostream>
#include <fstream>

GlobalArrayTransfer *gat;

namespace Verdandi
{
//! This class is the "model" interface for the simvascular/phasta flowsolver
//! to the Verdandi data assimilation library
//! It provides the necessary member functions for the Reduced Order Unscented Kalman Filter driver
class SimvascularVerdandiModel: public VerdandiBase
{
public:
	//! The numerical type (e.g., double).
	typedef double value_type;
	//! Type of the state error variance.
	typedef Matrix<double, General, PETScMPIDense> state_error_variance;
	//! Type of a row of the background error variance.
	typedef Vector<double, PETScPar> state_error_variance_row;
	//! Type of the state/observation crossed matrix.
	typedef Matrix<double, General, PETScMPIDense> matrix_state_observation;
	//! Type of the tangent linear model.
	typedef Matrix<double, General, PETScMPIDense> tangent_linear_operator;
	//! Type of the state vector.
	typedef Vector<double, PETScPar> state;

	//! Type of the reduced state error variance
	typedef Matrix<double, General, RowMajor, MallocAlloc<double> > reduced_state_error_variance;

protected:

	int Nparameter_;
	int Nstate_;
	int Nstate_local_;

	int rank_;
	int numProcs_;

	MPI_Comm iNewComm_C_;

	state duplicated_state_;

	int solveReturn_;

	//! Background error variance.
	double state_error_variance_value_;

	ofstream param_out_;

public:

	// Constructor and destructor.
	SimvascularVerdandiModel();
	~SimvascularVerdandiModel();
	void Initialize(string configuration_file);

	void Initialize();

	void InitializeFirstStep();
	void InitializeStep();
	void Finalize();

	// Processing.
	void Forward();
	void ForwardFinalize();
	bool HasFinished() const;

	// Operators.
	void ApplyOperator(state& x,
			bool forward = false, bool preserve_state = true);

	// Access methods.
	double GetTime() const;
	void SetTime(double time);
	int GetNstate() const;
	int GetLocalNstate() const;
	state& GetState();
	void StateUpdated();

	int GetNumProcs() const;
	int GetRank() const;

	// Errors.
	template <class L_matrix, class U_matrix>
	void GetStateErrorVarianceSqrt(L_matrix& L, U_matrix& U);

	string GetName() const;
	void Message(string message);

};


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_SimvascularVerdandiModel_HXX
#endif
