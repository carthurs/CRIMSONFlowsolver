#ifndef SIMVASCULAROBSERVATIONMANAGER_HXX

#include <limits>
#include <iostream>
#include <fstream>

#include "common_c.h"
#include "cvSolverIO.h"
#include "SimvascularVerdandiModel.hxx"

//#include "CGAL/Simple_cartesian.h"
//#include "CGAL/AABB_tree.h"
//#include "CGAL/AABB_traits.h"
//#include "CGAL/Polyhedron_3.h"
//#include "CGAL/AABB_polyhedron_segment_primitive.h"
//
//typedef CGAL::Simple_cartesian<double> K;
//typedef K::Point_3 Point;
//typedef K::Plane_3 Plane;
//typedef K::Vector_3 CGAL_Vector;
//typedef K::Segment_3 Segment;
//typedef CGAL::Polyhedron_3<K> Polyhedron;
//typedef CGAL::AABB_polyhedron_segment_primitive<K,Polyhedron> Primitive;
//typedef CGAL::AABB_traits<K, Primitive> Traits;
//typedef CGAL::AABB_tree<Traits> Tree;
//typedef Tree::Object_and_primitive_id Object_and_primitive_id;
//typedef Tree::Primitive_id Primitive_id;

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
	//! Period with which available observations are actually loaded.
	int Nskip_;
	//! First time at which observations are available.
	double initial_time_;
	//! Final time at which observations are available.
	double final_time_;

	//! Do we have simulation results as synthetic data?
	int synthetic_data_;

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

    int* nodes_uniq_;


    Vector<int> StateObsIndex_;
    Vector<int> DataArraysObsIndex_;

	int rank_;
	int numProcs_;

	/*** Cross-sectional flow observation ***/
	vtkSmartPointer<vtkPoints> geom_points_;
	vtkSmartPointer<vtkIdList> geom_ids_;
	vtkSmartPointer<vtkUnstructuredGrid> geom_UGrid_;

	vtkSmartPointer<vtkPlane> geom_plane_;
	vtkSmartPointer<vtkCutter> geom_cutter_;
	vtkSmartPointer<vtkConnectivityFilter> geom_connectivity_;

	vtkSmartPointer<vtkDoubleArray> geom_vel_array_;

	ofstream flow_out_;

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

	void GetObservationFlow(observation& observation);

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

	template<class state>
	void ApplyOperatorFlow(const state& x) ;

	template<class state>
	void ObserveFlow(const state& x, observation& y) const;

	double GetErrorVariance(int i, int j) const;
	const error_variance& GetErrorVariance() const;
	const error_variance& GetErrorVarianceInverse() const;

	void loadrestart(int timeindex, double* soln, double* acc, double* disp);

	string GetName() const;
	void Message(string message);
};

} // namespace Verdandi.

class is_near
{
public:
	bool operator() (vector<double> first, vector<double> second)
	{
		return (fabs(first[0] - second[0]) < std::numeric_limits<double>::epsilon() &&
			   fabs(first[1] - second[1]) < std::numeric_limits<double>::epsilon() &&
			   fabs(first[2] - second[2]) < std::numeric_limits<double>::epsilon() );
//		return (sqrt((first[0]-second[0])*(first[0]-second[0]) +
//				(first[1]-second[1])*(first[1]-second[1]) +
//				(first[2]-second[2])*(first[2]-second[2])) < std::numeric_limits<double>::epsilon());
	}
};

bool compare_firstcoord (vector<double> &first, vector<double> &second)
{
	return first[0] < second[0];
}

bool compare_secondcoord (vector<double> &first, vector<double> &second)
{
	return first[1] < second[1];
}

bool compare_thirdcoord (vector<double> &first, vector<double> &second)
{
	return first[2] < second[2];
}

#define SIMVASCULAROBSERVATIONMANAGER_HXX
#endif
