#ifndef VERDANDI_FILE_MODEL_SimvascularVerdandiModel_HXX
#define VERDANDI_FILE_MODEL_SimvascularVerdandiModel_HXX

// this brings in the common block global variables
#include "common_c.h"

// includes for the PHASTA Fortran functions
#include "partition.h"
#include "input_fform.h"
#include "input.h"
#include "proces.h"
#include "itrdrv.h"
#include "estimation_helpers.h"

#include "SimvascularGlobalArrayTransfer.h"

#include "mpi.h"

#ifdef intel
#include <direct.h>
#else
#include <unistd.h>
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

//! This class organizes the pointers to the PHASTA arrays
//! into augmented state components
class SimvascularAugStatePart
{
public:

	//! Constructor.
	SimvascularAugStatePart();

	//! Destructor
	~SimvascularAugStatePart();

	//! Initialize with a name
	void Initialize(const std::string& setname);

	//! Adds a pointer to the array of pointers
	void addDataPointer(double* data_pointer);

	//! Adds to the array of boolean flags for the estimated variables
	void addIsEstimated(bool val);

	//! Sets a specific pointer in the array of pointers
	void setDataPointer(int position, double* data_pointer);

	//! Sets the data for a specific variable.
	void setData(int position, double val);

	//! Sets the estimated-variable-flag for a specific variable
	void setIsEstimated(int position, bool val);

	//! Sets the constant that the output of getData is multiplied with
	void setPremulConstant(double val);

	//! Returns the data pointer to a specific variable
	double * getDataPointer(int position);

	//! Returns the value of a specific variable
	double getData(int position);

	//! Returns the constant that the output of getData is multiplied with
	double getPremulConstant();

	//! Returns the estimated-variable-flag for a specific variable
	bool getIsEstimated(int position);

	//! Returns the number of variables
	std::size_t getSize();

	//! Returns the number of estimated variables
	unsigned int getNumEstimated();

	//! Returns the index of the first estimated variable
    int getFirstEstimated();

    //! Returns the name
	std::string getName();

	//! Empty all the arrays and sets the multiplying constant to 1.0
	void Clear();

protected:

	//! Vector of pointers to PHASTA Fortran data
	std::vector<double*> pointers_;

	//! Vector of flags denoting if quantity is estimated
	std::vector<bool> is_estimated_;

	//! Name of the augmented state component
	std::string name_;

	//! The constant that multiplies the
	double premulconst_;
};

namespace Verdandi
{
//! This class is the "model" interface for the simvascular/phasta flowsolver
//! to the Verdandi data assimilation library.
//! It provides the necessary member functions for the Reduced Order Unscented Kalman Filter driver.
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

	//! Organizes the augmented state variables in a way that is expected by Verdandi.
	state duplicated_state_;

	//! Vector of distributed augmented state components
	std::vector<SimvascularAugStatePart> dstrb_parts_;

	//! Vector of non-distributed augmented state components
	std::vector<SimvascularAugStatePart> shared_parts_;


	//! Estimation error covariance diagonal values
	std::vector<double> state_error_variance_value_;


	//! Flags to denote that individual Windkessel-RCR R2 are estimated
	std::vector<int> cp_rcr_include_resistance_;

	//! Flags to denote that individual Windkessel-RCR C are estimated
	std::vector<int> cp_rcr_include_compliance_;

	//! Flags to denote that individual Windkessel-RCR R1 are estimated
	std::vector<int> cp_rcr_include_prox_resistance_;


	//! Pointer to the singleton instance of SimvascularGlobalArrayTransfer
	SimvascularGlobalArrayTransfer *gat;


	//! Size of the estimated portion of duplicated_state_
	int Nreduced_;

	//! Size of the entire duplicated_state_
	int Nstate_;

	//! Size of the local part of duplicated_state_
	int Nstate_local_;

	//! Starting index in duplicated_state_ of the estimated portion
	int state_reduced_start_local_;

	//! Size of the non-distributed augmented state variables
	int shared_parts_size_;


	//! Rank of the MPI process
	int rank_;

	//! Number of MPI processes
	int numProcs_;


	//! Flag to denote that vessel wall parameters are estimated
	int nreduced_has_wall_parameters_;

	//! Flag to denote that LPN parameters are estimated
	int nreduced_has_coupled_parameters_;

	//! Flag to denote that Windkessel-RCR R2 are estimated
	int cp_rcr_estimate_resistance_;

	//! Flag to denote that Windkessel-RCR C are estimated
	int cp_rcr_estimate_compliance_;

	//! Flag to denote that Windkessel-RCR R1 are estimated
	int cp_rcr_estimate_prox_resistance_;

	//! Flag to denote that Windkessel-RCR ground pressure is estimated
	int cp_rcr_estimate_pout_;

	//! Flag to denote that Windkessel-RCR pressure state variables are estimated
	int cp_rcr_estimate_pstates_;

    //! MPI communicator
	MPI_Comm iNewComm_C_;


	//! File stream for outputing estimated parameters
	ofstream Eoutfile_;

	//! File stream for reading in estimated parameters
	ifstream Einfile_;

public:

	//! Constructor.
	SimvascularVerdandiModel();

	//! Destructor.
	~SimvascularVerdandiModel();

	//! Initializes the model.
	void Initialize(string configuration_file);

	//! Initializes the model.
	void Initialize();

	//! Organizes the array pointers from simvascular
	void BuildAugmentedState();

	//! Initializes the first time step for the model.
	void InitializeFirstStep();

	//! Initializes the current time step for the model.
	void InitializeStep();

	//! Finalizes the current time step for the model.
	void Finalize();

	// Processing.

	//! Advances one step forward in time.
	void Forward();

	//! Finalizes the current time step
	void FinalizeStep();

	//! Checks whether the model has finished.
	bool HasFinished() const;

	// Operators.

	//! Applies the model (forward dynamics) operator to a given state.
	void ApplyOperator(state& x, bool forward = false, bool preserve_state = true);

	// Access methods.

	//! Returns the current time.
	double GetTime() const;

	//! Sets the time of the model to a given time.
	void SetTime(double time);

	//! Returns the state vector size.
	int GetNstate() const;

	//! Returns the size of the local (on-processor) state vector.
	int GetLocalNstate() const;

	//! Returns the starting index of the local reduced state vector
	int GetLocalReducedStart() const;

	//! copies internal state to "duplicated state" and returns a reference to it
	state& GetState();

	//! Updates the internal state from the duplicated state
	void StateUpdated();

	//! Returns the number of MPI
	int GetNumProcs() const;

	//! Returns the name of the class.
	int GetRank() const;

	// Errors.

	//! Returns a decomposition of the initial estimation error covariance matrix
	template <class L_matrix, class U_matrix>
	void GetStateErrorVarianceSqrt(L_matrix& L, U_matrix& U);


	//! Write only the estimated parameters to file
	void WriteEstimates();

	//! Returns the name of the class.
	string GetName() const;

	//! Receives and handles a message.
	void Message(string message);

};


} // namespace Verdandi.

#endif
