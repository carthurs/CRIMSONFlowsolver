#ifndef VERDANDI_FILE_MODEL_SimvascularVerdandiModel_HXX

#include "common_c.h"

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
class SimvascularAugStatePart
{
protected:

	std::vector<double*> pointers_;
	std::vector<bool> is_estimated_;
	std::vector<double> state_error_variance_value_;

	std::string name_;

	double premulconst_;

public:

	SimvascularAugStatePart();
	~SimvascularAugStatePart();

	void Initialize(const std::string& setname);

	void addDataPointer(double* data_pointer);
	void addIsEstimated(bool val);

	void setDataPointer(int position, double* data_pointer);
	void setData(int position, double val);
	void setIsEstimated(int position, bool val);
	void setPremulConstant(double val);

	double * getDataPointer(int position);
	double getData(int position);

	double getPremulConstant();

	bool getIsEstimated(int position);
	std::size_t getSize();
	std::size_t getNumEstimated();
    int getFirstEstimated();

	std::string getName();

	void Clear();

};

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

	int Nreduced_;
	int Nstate_;
	int Nstate_local_;
	int state_reduced_start_local_;
	int shared_parts_size_;

	int rank_;
	int numProcs_;

	int nreduced_has_wall_parameters_;
	int nreduced_has_coupled_parameters_;

	int cp_rcr_estimate_resistance_;
	int cp_rcr_estimate_compliance_;
	int cp_rcr_estimate_prox_resistance_;
	int cp_rcr_estimate_pout_;
	int cp_rcr_estimate_pstates_;

	std::vector<int> cp_rcr_include_resistance_;
	std::vector<int> cp_rcr_include_compliance_;
	std::vector<int> cp_rcr_include_prox_resistance_;

	std::vector<int> cp_rcr_face_grouping_;

	MPI_Comm iNewComm_C_;

	std::vector<double> state_error_variance_value_;

	SimvascularGlobalArrayTransfer *gat;

	state duplicated_state_;

	std::vector<SimvascularAugStatePart> dstrb_parts_;
	std::vector<SimvascularAugStatePart> shared_parts_;

	ofstream Eoutfile_;
	ifstream Einfile_;

public:

	// Constructor and destructor.
	SimvascularVerdandiModel();
	~SimvascularVerdandiModel();

	void Initialize(string configuration_file);
	void Initialize();
	void BuildAugmentedState();

	void InitializeFirstStep();
	void InitializeStep();
	void Finalize();

	// Processing.
	void Forward();
	void FinalizeStep();
	bool HasFinished() const;

	// Operators.
	void ApplyOperator(state& x, bool forward = false, bool preserve_state = true);

	// Access methods.
	double GetTime() const;
	void SetTime(double time);
	int GetNstate() const;
	int GetLocalNstate() const;
	int GetLocalReducedStart() const;
	state& GetState();
	void StateUpdated();

	int GetNumProcs() const;
	int GetRank() const;

	// Errors.
	template <class L_matrix, class U_matrix>
	void GetStateErrorVarianceSqrt(L_matrix& L, U_matrix& U);

	// Output
	void WriteEstimates();

	string GetName() const;
	void Message(string message);

};


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_SimvascularVerdandiModel_HXX
#endif
