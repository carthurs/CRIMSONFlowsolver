#ifndef VERDANDI_FILE_MODEL_SimvascularVerdandiModel_HXX

#include "partition.h"
#include "phSolver.h"
#include "SCField.h"

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
	typedef Matrix<double> state_error_variance;
	//! Type of a row of the background error variance.
	typedef Vector<double> state_error_variance_row; // not needed for ROUKF but required for compile
	//! Type of the state/observation crossed matrix.
	typedef Matrix<double> matrix_state_observation;
	//! Type of the tangent linear model.
	typedef Matrix<double> tangent_linear_operator;  // not needed for ROUKF but required for compile
	//! Type of the state vector.
	typedef Vector<double> state;

	int Nparameter_;
	int Nstate_;
	int Nobservation_;

	state duplicated_state_;

	phSolver* phS;
	int solveReturn;


//public:

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
	int GetNfull_state() const;
	void GetState(state& state);
	void SetState(state& state);
	void GetFullState(state& state);
	void SetFullState(state& state);

	// Errors.
	void GetStateErrorVarianceSqrt(state_error_variance& L,
			state_error_variance& U);

	string GetName() const;
	void Message(string message);
};


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_SimvascularVerdandiModel_HXX
#endif
