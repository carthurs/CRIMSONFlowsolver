#ifndef SIMVASCULAROBSERVATIONMANAGER_HXX
#define SIMVASCULAROBSERVATIONMANAGER_HXX

#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "mpi.h"

// VTK includes
#define VTK_EXCLUDE_STRSTREAM_HEADERS
#include "vtkIdList.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkPolyDataWriter.h"
#include "vtkSmartPointer.h"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkConnectivityFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkDataSetSurfaceFilter.h"

// includes for the PHASTA Fortran functions
#include "distmeas.h"
#include "elmdist.h"

namespace Verdandi {

//! This is the base class for the different types of observations
class SimvascularObservationType {

public:

	//! Type of the observation vector.
    typedef Vector<double> observation;

    //! Type of the state vector.
    typedef Vector<double, PETScPar> state;

protected:

	//! Indices for nodal observation operators
	std::vector<int> state_obs_index_;

	//! vector of measurement error variances (diagonal entries)
	std::vector<double> error_variance_value_;

	//! name of the type of observation
	std::string name_;

	//! File name for the data file
	std::string data_file_name_;

	//! Name of directory with data file
	std::string data_directory_;

	//! Name of file to which observations will be saved
	std::string saved_obs_file_name_;

	//! Input file stream
	std::ifstream data_file_in_;

	//! Output file stream
	std::ofstream data_file_out_;

	//! Pointer to the single instance of SimvascularGlobalArrayTransfer
	SimvascularGlobalArrayTransfer *gat;

	//! Number of observations on this process
	int Nobservation_local_;

	//! Rank of MPI Process
	int rank_;

	//! Number of MPI Processes
	int numProcs_;

	//! MPI communicator
	MPI_Comm iNewComm_C_;

	//! Flag for loading data from a .dat file
	bool load_data_from_datfile_;

	//! Flag to denote whether or not the observations are distributed among processes
	bool isDistributed_;

public:

	//! Constructor
	SimvascularObservationType();

	//! Destructor
	~SimvascularObservationType();

	//! Initialization
	virtual void Initialize(std::string name, const SimvascularAugStatePart &aug_state_part, VerdandiOps &configuration) = 0;

	//! Apply the observation operator
	virtual void ApplyOperator(const state& x, observation& Hx1, observation& Hx2, int obs_start_index) const = 0;

	//! Load the measurements from file
	void LoadData(int linetoread, Vector<double>& dataarray, int obs_start_index);

	//! Save the observations to data
	void SaveObservations(const state& x);

	//! Get the number of local observations
	int getNobservation_local() const;

	//! Get a measurement error variance (diagonal entry)
	double getErrorVarianceValue(int ind) const;

	//! Get the name of the observation type
	std::string getName() const;

};

class SimvascularNodalSolutionObservation : public SimvascularObservationType {
protected:
public:

	//! Constructor
	SimvascularNodalSolutionObservation();

	//! Destructor
	~SimvascularNodalSolutionObservation();

	//! Initialization
	void Initialize(std::string name, const SimvascularAugStatePart &aug_state_part, VerdandiOps &configuration);

	//! Apply the observation operator
	void ApplyOperator(const state& x, observation& Hx1, observation& Hx2, int obs_start_index) const;

};

class SimvascularNodalDisplacementObservation : public SimvascularObservationType {
protected:
public:

	//! Constructor
	SimvascularNodalDisplacementObservation();

	//! Destructor
	~SimvascularNodalDisplacementObservation();

	//! Initialization
	void Initialize(std::string name, const SimvascularAugStatePart &aug_state_part, VerdandiOps &configuration);

	//! Apply the observation operator
	void ApplyOperator(const state&x, observation& Hx1, observation& Hx2, int obs_start_index) const;

};

class SimvascularDistanceObservation : public SimvascularObservationType {
protected:

	//! Indices in the duplicated state vector for the distance field
	std::vector<int> disp_duplicated_state_index_;
public:

	//! Constructor
	SimvascularDistanceObservation();

	//! Destructor
	~SimvascularDistanceObservation();

	//! Initialization
	void Initialize(std::string name, const SimvascularAugStatePart &aug_state_part, VerdandiOps &configuration);

	//! Apply the observation operator
	void ApplyOperator(const state& x, observation& Hx1, observation& Hx2, int obs_start_index) const;

};

class SimvascularFlowPressObservation : public SimvascularObservationType {
protected:

	//! Origins of the cut planes
	std::vector<Seldon::Vector<double> > csobs_origins_;

	//! Normals of the cut planes
	std::vector<Seldon::Vector<double> > csobs_normals_;

	//! Radii associated with the cut planes
	std::vector<double> csobs_radii_;

    //! Indices in the duplicated state for the solution field
	std::vector<int> sol_duplicated_state_index_;

	// pointers to VTK objects for reslicing and resampling solution field
	vtkSmartPointer<vtkPoints> geom_points_;
	vtkSmartPointer<vtkPoints> geom_points_def_;
	vtkSmartPointer<vtkIdList> geom_ids_;
	vtkSmartPointer<vtkUnstructuredGrid> geom_UGrid_;
	vtkSmartPointer<vtkUnstructuredGrid> geom_UGrid_def_;
	vtkSmartPointer<vtkDataSetSurfaceFilter> geom_surface_def_;
	vtkSmartPointer<vtkPolyDataWriter> geom_writer_;
	vtkSmartPointer<vtkPlane> geom_plane_;
	vtkSmartPointer<vtkCutter> geom_cutter_;
	vtkSmartPointer<vtkCutter> geom_cutter_alt_;
	vtkSmartPointer<vtkPolyDataWriter> geom_writers_;

	//! vector of pointers to VTK plane objects representing cross-section cuts
	std::vector <vtkSmartPointer<vtkPlane> > geom_planes_;

	//! vector of pointers to VTK cutter objects for the cross-section cuts
	std::vector <vtkSmartPointer<vtkCutter> > geom_cutters_;

	//! vector of vectors storing distances from resampled cells to origin points
	std::vector<std::vector<double> > distances_fromorigin_;

	//! number of flow observations
    int Nobservation_flow_;

    //! number of average pressure observations
    int Nobservation_avgpressure_;

    //! number of area observations
    int Nobservation_area_;

public:

    //! Constructor
	SimvascularFlowPressObservation();

	//! Destructor
	~SimvascularFlowPressObservation();

	//! Initialization
	void Initialize(std::string name, const SimvascularAugStatePart &aug_state_part, VerdandiOps &configuration);

	//! Apply the observation operator
	void ApplyOperator(const state& x, observation& Hx1, observation& Hx2, int obs_start_index) const;

};

///////////////////////////////////
// SIMVASCULAROBSERVATIONMANAGER //
///////////////////////////////////

//! This class is an observation manager for the Simvascular flowsolver.
class SimvascularObservationManager: public VerdandiBase {

public:

	//! Type of the tangent linear operator
	typedef Matrix<double, General, PETScMPIAIJ> tangent_linear_operator;

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

	std::vector <SimvascularObservationType*> observations_dstrb_;
	std::vector <SimvascularObservationType*> observations_single_;

	//! Directory that stores the forward results.
	std::string data_directory_;

	//! Total number of observations at current time.
	int Nobservation_;

	//! Local number of observations at current time.
	int Nobservation_local_;

	//! Number of cross-sectional flow observations
	int Nobservation_flow_;

	//! Number of cross-sectional average pressure observations
	int Nobservation_avgpressure_;

	//! Period with which available observations are actually loaded.
	int Nskip_;

	//! First time at which observations are available.
	double initial_time_;

	//! Final time at which observations are available.
	double final_time_;

	//! Flag that denotes whether we use simulation restarts as data
	int use_restarts_;

	//! Flag to carry out single node observations
    int execute_nodal_observations_;

    //! Flag to carry out distance observations
    int execute_distance_observations_;

	/*** Observation times ***/

	//! Requested time.
	double time_;

	//! Period length of data, if it is assumed to be periodic
	double data_period_;

	//! Flag to discard observation
	bool discard_observation_;

	 /*** Observation operator ***/

	//! Tangent operator matrix (H).
	tangent_linear_operator tangent_operator_matrix_;

	//! Observation error covariance matrix (R).
	error_variance error_variance_;

#ifdef VERDANDI_ROUKF_PARALLEL_INNOVATION
	//! Inverse of the diagonal observation error covariance matrix
	Vector<double, PETScPar> error_variance_inverse_diag_;
#endif

	//! Inverse of the observation error covariance matrix
	error_variance error_variance_inverse_;

    /*** Model domain ***/

    //! The size of a model state.
    int Nstate_model_;

    //! Rank of MPI Process
	int rank_;

	//! Number of MPI Processes
	int numProcs_;

    /*** Observation data (measurements) ***/

    //! Lower value for linear interpolation
    int current_lower_bound_;

    //! Upper value for linear interpolation
    int current_upper_bound_;

    //! pointer to the single instance of SimvascularGlobalArrayTransfer
	SimvascularGlobalArrayTransfer *gat;

    //! Arrays for linear interpolation
    Vector<double> dataarrays_lower_;

    //! Arrays for linear interpolation
    Vector<double> dataarrays_upper_;


	/*** MPI ***/

	//! MPI communicator
	MPI_Comm iNewComm_C_;

public:

	//! Default constructor.
	SimvascularObservationManager();

	//! Destructor.
	~SimvascularObservationManager();

	/*** Initialization ***/

	//! Initializes the observation manager.
	template<class Model>
	void Initialize(const Model& model, string configuration_file);

	//! Sets the time of observations to be loaded.
	template<class Model>
	void SetTime(const Model& model, double time);

    /*** Methods for observation data ***/

	//! Indicates if some observations are available at a given time.
	bool HasObservation() const;

	//! Indicates if some observations are available at current time.
	bool HasObservation(double time);

	//! Activates or deactivates the option 'discard_observation'.
	void DiscardObservation(bool discard_observation);

	//! Returns the number of available observations.
	int GetNobservation() const;

	//! Returns the size of the local (on-processor) state vector.
	int GetLocalNobservation() const;

	//! Returns the observations. Not implemented, as we prefer GetInnovation
	void GetObservation(observation& observation);
	//void loadrestart(int timeindex, double* soln, double* acc, double* disp);

	//! Load data
	void LoadObservationSingleLocal(int timeindex, Vector<double>& dataarray);

	//! Saves observations in file
	template<class state>
	void SaveObservationSingleLocal(const state& x);




	/*** Operators ***/

	//! Applies the observation operator to a given vector.
	template<class state>
	void ApplyOperator(const state& x, observation& y) const;

	//! Applies the observation operator to a given vector (local to process).
	template<class state>
	void ApplyOperatorLocal(const state& x, observation& Hx1, observation& Hx2);

	//! Applies the flow/ avg. pressure observation operator to a given state.
	template<class state>
	void ApplyOperatorFlow(const state& x, observation& Hx);

	//! Gets the innovation.
	template<class state>
	void GetInnovation(const state& x, observation& innovation);

#if defined(VERDANDI_ROUKF_PARALLEL_INNOVATION)
	//! Gets the innovation.
	template<class state>
	void GetInnovation(const state& x,state& innovation_p_orig, state& innovation_p_fe);
#else
	template<class state>
	void GetInnovation(const state& x,observation& innovation_orig, observation& innovation_fe);
#endif

	//! Return an observation error variance.
	double GetErrorVariance(int i, int j) const;

	//! Returns the observation error variance.
	const error_variance& GetErrorVariance() const;

	//! Returns the inverse of the error variance matrix.
	const error_variance& GetErrorVarianceInverse() const;

	//! Returns the name of the class.
	string GetName() const;

	void Message(string message);
};

} // namespace Verdandi.

#endif
