#ifndef SIMVASCULAROBSERVATIONMANAGER_HXX

namespace Verdandi {

///////////////////////////////////
// SIMVASCULAROBSERVATIONMANAGER //
///////////////////////////////////

//! This class is a template of observation manager.
class SimvascularObservationManager: public VerdandiBase {
public:
#ifdef VERDANDI_TANGENT_LINEAR_OPERATOR_SPARSE
	//! Type of the tangent linear operator.
	typedef Matrix<double, General, RowSparse> tangent_linear_operator;
#else
	//! Type of the tangent linear operator.
	typedef Matrix<double> tangent_linear_operator;
#endif
	//! Type of a row of the tangent linear operator.
	typedef Vector<double> tangent_linear_operator_row;

#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
	//! Type of the observation error covariance matrix.
	typedef Matrix<double, General, RowSparse> error_variance;
#else
	//! Type of the observation error covariance matrix.
	typedef Matrix<double> error_variance;
#endif

	//! Type of the observation vector.
	typedef Vector<double> observation;

protected:

	/*** Observation file structure ***/

	//! Directory that stores the forward results.
	string data_directory_;


	//! Total number of observations at current time.
	int Nobservation_;
	//! Period with which available observations are actually loaded.
	int Nskip_;
	//! First time at which observations are available.
	double initial_time_;
	//! Final time at which observations are available.
	double final_time_;

	/*** Observation times ***/

	//! Requested time.
	double time_;

	bool discard_observation_;

	 /*** Observation operator ***/

	//! Tangent operator matrix (H).
	tangent_linear_operator tangent_operator_matrix_;

	//! Observation error variance.
	double error_variance_value_;
	//! Observation error covariance matrix (R).
	error_variance error_variance_;
	//! Inverse of the observation error covariance matrix (R).
	error_variance error_variance_inverse_;

    /*** Model domain ***/

    //! The size of a model state.
    int Nstate_model_;

    int current_upper_bound_;
    int current_lower_bound_;

    double* dataarrays_lower_;
    double* dataarrays_upper_;
    double* soln_lower_;
    double* soln_upper_;
    double* acc_lower_;
    double* acc_upper_;
    double* disp_lower_;
    double* disp_upper_;

    int isize_solution_;
    int isize_displacement_;
    int isize_nshg_;
    int isize_nshguniq_;

    int* linobs_soln_;
    int* linobs_acc_;
    int* linobs_disp_;

    int* nodes_uniq_;

    Vector<int> StateObsIndex_;
    Vector<int> DataArraysObsIndex_;

public:
	// Constructors and destructor.
	SimvascularObservationManager();
	~SimvascularObservationManager();

	// Initialization.
	template<class Model>
	void Initialize(const Model& model, string configuration_file);

	template<class Model>
	void SetTime(const Model& model, double time);

	void DiscardObservation(bool discard_observation);

	/////////////////
	// OBSERVATION //
	/////////////////

	void GetObservation(observation& observation);

	////////////////
	// INNOVATION //
	////////////////

	template<class state>
	void GetInnovation(const state& x, observation& innovation);

	////////////
	// ACCESS //
	////////////

	bool HasObservation() const;
	bool HasObservation(double time);
	int GetNobservation() const;

	///////////////
	// OPERATORS //
	///////////////

	template<class state>
	void ApplyOperator(const state& x, observation& y) const;

	double GetErrorVariance(int i, int j) const;
	const error_variance& GetErrorVariance() const;
	const error_variance& GetErrorVarianceInverse() const;

	void loadrestart(int timeindex, double* soln, double* acc, double* disp);

	string GetName() const;
	void Message(string message);
};

} // namespace Verdandi.

#define SIMVASCULAROBSERVATIONMANAGER_HXX
#endif
