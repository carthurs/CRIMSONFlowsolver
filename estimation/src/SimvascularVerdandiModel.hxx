#ifndef VERDANDI_FILE_MODEL_SimvascularVerdandiModel_HXX

#include "partition.h"
#include "phSolver.h"
#include "SCField.h"
#include "estimation_helpers.h"
#include "common_c.h"

#include <iostream>
#include <fstream>

phSolver* phS;

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

protected:

	int Nparameter_;
	int Nstate_;
	int Nstate_local_;

	int rank_;
	int numProcs_;

	state duplicated_state_;

	int solveReturn;

	const SCField *node_field_,*soln_field_,*acc_field_,*disp_field_;

	//! Background error variance.
	double state_error_variance_value_;

	ofstream param_out_;

public:

	// Constructor and destructor.
	SimvascularVerdandiModel();
	~SimvascularVerdandiModel();
	void Initialize(string configuration_file);

	void Initialize(int argc, char * argv[]);

	void InitializeFirstStep();
	void InitializeStep();

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

	// Errors.
	template <class L_matrix, class U_matrix>
	void GetStateErrorVarianceSqrt(L_matrix& L, U_matrix& U);

	string GetName() const;
	void Message(string message);

};


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_SimvascularVerdandiModel_HXX
#endif
