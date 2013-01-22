#ifndef SIMVASCULAROBSERVATIONMANAGER_HXX

#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#include "common_c.h"
#include "cvSolverIO.h"
#include "distmeas.h"

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
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"

namespace Verdandi {

///////////////////////////////////
// SIMVASCULAROBSERVATIONMANAGER //
///////////////////////////////////

//! This class is an observation manager for the Simvascular flowsolver.
class SimvascularObservationManager: public VerdandiBase {

public:

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

    int ignore_nodal_observations_;
    int ignore_distance_observations_;

	/*** Observation times ***/

	//! Requested time.
	double time_;

	double data_period_;

	bool discard_observation_;

	 /*** Observation operator ***/

	//! Tangent operator matrix (H).
	tangent_linear_operator tangent_operator_matrix_;
	//! Indices for simple observation operator
	Vector<int> StateObsIndex_;
	Vector<int> DataArraysObsIndex_;

	//! Observation error variance.
	double error_variance_value_;
	//! Observation error covariance matrix (R).
	error_variance error_variance_;
	//! Inverse of the observation error covariance matrix (R).
	error_variance error_variance_inverse_;

    /*** Model domain ***/

    //! The size of a model state.
    int Nstate_model_;

    /*** Observation data (measurements) ***/
    int current_upper_bound_;
    int current_lower_bound_;

    //! Arrays for linear interpolation
    double* dataarrays_lower_;
    double* dataarrays_upper_;
    double* soln_lower_;
    double* soln_upper_;
    double* acc_lower_;
    double* acc_upper_;
    double* disp_lower_;
    double* disp_upper_;

    //! Size of internal arrays
    int isize_solution_;
    int isize_displacement_;
    //! Number of global shape functions
    int isize_nshg_;
    //! Number of global shape functions only in master images
    int isize_nshguniq_;

	/*** Cross-sectional flow and pressure observation ***/
	vector<Seldon::Vector<double> > csobs_origins_;
	vector<Seldon::Vector<double> > csobs_normals_;
	vector<double> csobs_radii_;

	vtkSmartPointer<vtkPoints> geom_points_;
	vtkSmartPointer<vtkIdList> geom_ids_;
	vtkSmartPointer<vtkUnstructuredGrid> geom_UGrid_;

	vtkSmartPointer<vtkPlane> geom_plane_;
	vtkSmartPointer<vtkCutter> geom_cutter_;

	vector <vtkSmartPointer<vtkPlane> > geom_planes_;
	vector <vtkSmartPointer<vtkCutter> > geom_cutters_;

	vector<vector<double> > distances_fromorigin_;

	/***  ***/

	/*** File handling ***/
    string obsfilename_part_;
    string obsfilename_single_;

	ifstream obs_in_part_;
	ifstream obs_in_single_;

	ofstream obs_out_part_;
	ofstream obs_out_single_;
	//ofstream flow_out_;

	/*** MPI ***/
	int rank_;
	int numProcs_;

	int *obs_recvcount_;

public:

	SimvascularObservationManager();
	~SimvascularObservationManager();

	/*** Initialization ***/
	template<class Model>
	void Initialize(const Model& model, string configuration_file);
	template<class Model>
	void SetTime(const Model& model, double time);
	void InitializeFiles();

    /*** Methods for observation data ***/
	bool HasObservation() const;
	bool HasObservation(double time);
	void DiscardObservation(bool discard_observation);
	int GetNobservation() const;
	void GetObservation(observation& observation);
	void loadrestart(int timeindex, double* soln, double* acc, double* disp);
	void LoadObservationSingleLocal(int timeindex, double* dataarray);
	template<class state>
	void SaveObservationSingleLocal(const state& x);

	/*** Operators ***/
	template<class state>
	void ApplyOperator(const state& x, observation& y) const;
	template<class state>
	void ApplyOperatorLocal(const state& x, observation& Hx);
	template<class state>
	void ApplyOperatorFlow(const state& x, observation& Hx);
	template<class state>
	void GetInnovation(const state& x, observation& innovation);
	double GetErrorVariance(int i, int j) const;
	const error_variance& GetErrorVariance() const;
	const error_variance& GetErrorVarianceInverse() const;

	string GetName() const;
	void Message(string message);
};

} // namespace Verdandi.

#define SIMVASCULAROBSERVATIONMANAGER_HXX
#endif
