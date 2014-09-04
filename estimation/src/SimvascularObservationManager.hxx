#ifndef SIMVASCULAROBSERVATIONMANAGER_HXX
#define SIMVASCULAROBSERVATIONMANAGER_HXX

#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

// this brings in the common block global variables
#include "common_c.h"

// includes for the PHASTA Fortran functions
#include "cvSolverIO.h"
#include "distmeas.h"
#include "elmdist.h"

#include "SimvascularGlobalArrayTransfer.h"

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

//class SimvascularObservationType {
//
//};

namespace Verdandi {

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

	//! Directory that stores the forward results.
	string data_directory_;

	//! Total number of observations at current time.
	int Nobservation_;

	//! Local number of observations at current time.
	int Nobservation_local_;

	//! Number of single node observations
	int Nobservation_nodal_;

	//! Number of distance observations
	int Nobservation_dist_;

	//! Number of area observations
	int Nobservation_area_;

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

	//! Indices for nodal observation operators
	Vector<int> StateObsIndex_;

	//! Indices for the data array
	Vector<int> DataArraysObsIndex_;

	//! Observation error variance.
	double error_variance_value_;

	//! diagonal error covariance value for nodal observations
	double error_variance_value_nodal_;

	//! diagonal error covariance value for distance observations
	double error_variance_value_dist_;

	//! diagonal error covariance value for pressure observations
	double error_variance_value_avgpress_;

	//! diagonal error covariance value for flow observations
	double error_variance_value_flow_;


	//! Observation error covariance matrix (R).
	error_variance error_variance_;

#ifdef VERDANDI_ROUKF_PARALLEL_INNOVATION
	//! Inverse of the diagonal observation error covariance matrix
	Vector<double, PETScPar> error_variance_inverse_diag_;
#endif
	//! Inverse of the observation error covariance matrix
	error_variance error_variance_inverse_;


    /*** Model domain ***/

	//! pointer to the single instance of SimvascularGlobalArrayTransfer
	SimvascularGlobalArrayTransfer *gat;

    //! The size of a model state.
    int Nstate_model_;

    /*** Observation data (measurements) ***/

    //! Upper value for linear interpolation
    int current_upper_bound_;

    //! Lower value for linear interpolation
    int current_lower_bound_;

    //! Arrays for linear interpolation
    Vector<double> dataarrays_lower_;

    //! Arrays for linear interpolation
    Vector<double> dataarrays_upper_;

    int isize_solution_;
    int isize_displacement_;

    //! Number of global shape functions
    int isize_nshg_;
    //! Number of global shape functions only in master images
    int isize_nshguniq_;

	/*** Cross-sectional flow and pressure observation ***/

    //! Origins of the cut planes
	vector<Seldon::Vector<double> > csobs_origins_;

	//! Normals of the cut planes
	vector<Seldon::Vector<double> > csobs_normals_;

	//! Radii associated with the cut planes
	vector<double> csobs_radii_;


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
	//vtkSmartPointer<vtkConnectivityFilter> geom_connectivity_;

	vector <vtkSmartPointer<vtkPlane> > geom_planes_;
	vector <vtkSmartPointer<vtkCutter> > geom_cutters_;
    //vector <vtkSmartPointer<vtkConnectivityFilter> > geom_connec_filters_;

	vector<vector<double> > distances_fromorigin_;

	/*** File handling ***/
    string obsfilename_part_;
    string obsfilename_single_;

	ifstream obs_in_part_;
	ifstream obs_in_single_;

	ofstream obs_out_part_;
	ofstream obs_out_single_;

	/*** MPI ***/

	//! MPI communicator
	MPI_Comm iNewComm_C_;

	//! Rank of MPI Process
	int rank_;

	//! Number of MPI Processes
	int numProcs_;

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
