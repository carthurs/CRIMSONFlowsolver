#ifndef SIMVASCULAROBSERVATIONMANAGER_CXX

#include "SimvascularObservationManager.hxx"

#include "cvSolverIO.h"

#include "common_c.h"

namespace Verdandi {

/////////////////////////////////
// CONSTRUCTORS AND DESTRUCTOR //
/////////////////////////////////

//! Default constructor.
SimvascularObservationManager::SimvascularObservationManager() {
    current_lower_bound_ = -1;
    current_upper_bound_ = -1;
}

//! Destructor.
SimvascularObservationManager::~SimvascularObservationManager() {
	// Operations to be performed when the object is destroyed.
	delete [] linobs_soln_;
	delete [] linobs_acc_;
	delete [] linobs_disp_;
}

////////////////////
// INITIALIZATION //
////////////////////

//! Initializes the observation manager.
/*!
 \param[in] model model.
 \param[in] configuration_file configuration file.
 \tparam Model the model type; e.g. ShallowWater<double>
 */
template<class Model>
void SimvascularObservationManager::Initialize(const Model& model,
		string configuration_file) {
	VerdandiOps configuration(configuration_file);

	int igeombc; /* file handle for geombc */
	int irestart; /* file handle for restart */
	int iarray[10];
	int ione = 1;
	int itwo = 2;
	int ithree = 3;
	char iformat[80];
	char filename[255];
	char restart_filename[255];
	char geombc_filename[255];

	strcpy(iformat, outpar.iotype);

	Nstate_model_ = model.GetNstate();

	configuration.SetPrefix("observation.");
	configuration.Set("data_directory", data_directory_);
	configuration.Set("Nskip", "v > 0", Nskip_);
	configuration.Set("initial_time", "", 0., initial_time_);
	configuration.Set("final_time", "", numeric_limits<double>::max(),
			final_time_);
	configuration.Set("error.variance", "v > 0", error_variance_value_);

	//
	// read the simple linear observation operator from the geombc files if found
	// these are distributed arrays
	//
	//
	// note that this geombc.dat.x should be the one in the current working folder
	//
	openfile_("geombc.dat.1", "read", &igeombc);

	// -------
	readheader_(&igeombc, "observation function solution?", (void*) iarray, &itwo, "integer",
			iformat);

	isize_solution_ = iarray[0] * iarray[1];
	isize_nshg_ = iarray[0];

	//cout << isize_solution_ << endl;

	linobs_soln_ = new int[isize_solution_];

	readdatablock_(&igeombc, "observation function solution?", (void*) linobs_soln_, &isize_solution_,
			"integer", iformat);

//	cout << "observation function solution:" << endl;
//
//	for (int jj = 0; jj < iarray[0]; jj++) {
//		for (int kk = 0; kk < iarray[1]; kk++) {
//			cout << linobs_soln_[kk*iarray[0]+jj] << " ";
//		}
//		cout << endl;
//	}
	// -------

	// -------
	readheader_(&igeombc, "observation function time derivative of solution?", (void*) iarray, &itwo, "integer",
				iformat);

	isize_solution_ = iarray[0] * iarray[1];

	linobs_acc_ = new int[isize_solution_];

	readdatablock_(&igeombc, "observation function time derivative of solution?", (void*) linobs_acc_, &isize_solution_,
			"integer", iformat);

//	cout << "observation function time derivative of solution:" << endl;
//
//	for (int jj = 0; jj < iarray[0]; jj++) {
//		for (int kk = 0; kk < iarray[1]; kk++) {
//			cout << linobs_acc_[kk*iarray[0]+jj] << " ";
//		}
//		cout << endl;
//	}
	// -------

	// -------
	readheader_(&igeombc, "observation function displacement?", (void*) iarray, &itwo, "integer",
			iformat);

	isize_displacement_ = iarray[0] * iarray[1];

	linobs_disp_ = new int[isize_displacement_];

	readdatablock_(&igeombc, "observation function displacement?", (void*) linobs_disp_, &isize_displacement_,
			"integer", iformat);

//	cout << "observation function displacement:" << endl;
//
//	for (int jj = 0; jj < iarray[0]; jj++) {
//		for (int kk = 0; kk < iarray[1]; kk++) {
//			cout << linobs_disp_[kk*iarray[0]+jj] << " ";
//		}
//		cout << endl;
//	}
	// -------


	//
	// read in the unique list of nodes
	//
    if (conpar.numnp > 1) {
    	readheader_(&igeombc, "local index of unique nodes?", (void*) iarray, &ione, "integer",
    				iformat);

    	isize_nshguniq_ = iarray[0];

    	nodes_uniq_ = new int[isize_nshguniq_];

    	readdatablock_(&igeombc, "local index of unique nodes?", (void*) nodes_uniq_, &isize_nshguniq_,
    				"integer", iformat);
    }

    //
    // allocate space for the data (right now it is the size of the full state)
    // these are distributed arrays
    // isize_solution_ and isize_displacement_ should be the same as the array sizes in restart,
    // or else we have had a problem
    //
    dataarrays_lower_ = new double[2*isize_solution_ + isize_displacement_];
    dataarrays_upper_ = new double[2*isize_solution_ + isize_displacement_];

    soln_lower_ = &dataarrays_lower_[0];
    soln_upper_ = &dataarrays_upper_[0];

    acc_lower_ = &dataarrays_lower_[isize_solution_];
    acc_upper_ = &dataarrays_upper_[isize_solution_];

    disp_lower_ = &dataarrays_lower_[2*isize_solution_];
    disp_upper_ = &dataarrays_upper_[2*isize_solution_];


//    soln_lower_ = new double[isize_solution_];
//    soln_upper_ = new double[isize_solution_];
//    acc_lower_ = new double[isize_solution_];
//    acc_upper_ = new double[isize_solution_];
//    disp_lower_ = new double[isize_displacement_];
//    disp_upper_ = new double[isize_displacement_];

    //
    // find out the number of observations
    // and do some preprocessing on the simple observation operator
    //
    int obsCounter = 0;

    for (int jj = 0; jj < isize_nshg_; jj++) {
    	for (int kk = 0; kk < 4; kk++) {
    		if (linobs_soln_[kk*isize_nshg_+jj] > 0) {
    			StateObsIndex_.PushBack(kk+4*jj);
    			DataArraysObsIndex_.PushBack(kk*isize_nshg_+jj);
    			obsCounter++;
    		}
    	}
    }

    for (int jj = 0; jj < isize_nshg_; jj++) {
    	for (int kk = 0; kk < 4; kk++) {
    		if (linobs_acc_[kk*isize_nshg_+jj] > 0) {
    			StateObsIndex_.PushBack(kk+4*jj + isize_nshg_*4);
    			DataArraysObsIndex_.PushBack(kk*isize_nshg_+jj + isize_solution_);
    			obsCounter++;
    		}
    	}
    }

    for (int jj = 0; jj < isize_nshg_; jj++) {
    	for (int kk = 0; kk < 3; kk++) {
    		if (linobs_disp_[kk*isize_nshg_+jj] > 0) {
    			StateObsIndex_.PushBack(kk+3*jj + isize_nshg_*4 + isize_nshg_*4);
    			DataArraysObsIndex_.PushBack(kk*isize_nshg_+jj + isize_solution_*2);
    			obsCounter++;
    		}
    	}
    }

    Nobservation_ = obsCounter;


    //
    // initial observation error variance matrix
    //

#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
    build_diagonal_sparse_matrix(Nobservation_, error_variance_value_,
    		error_variance_);
    build_diagonal_sparse_matrix(Nobservation_,
    		T(T(1) / error_variance_value_),
    		error_variance_inverse_);
#else
	error_variance_.Reallocate(Nobservation_, Nobservation_);
	error_variance_.SetIdentity();
	Mlt(error_variance_value_, error_variance_);
	error_variance_inverse_.Reallocate(Nobservation_, Nobservation_);
	error_variance_inverse_.SetIdentity();
	Mlt(double(double(1)/ error_variance_value_), error_variance_inverse_);
#endif

}

//! Sets the time of observations to be loaded.
/*!
 \param[in] model the model.
 \param[in] time a given time.
 */
template<class Model>
void SimvascularObservationManager::SetTime(const Model& model, double time) {
	if (time_ == time)
		return;

	time_ = time;
}

/////////////////
// OBSERVATION //
/////////////////

//! Activates or deactivates the option 'discard_observation'.
/*!
      \param[in] discard_observation if set to true, each observation will be
      used at most one time.
 */
void SimvascularObservationManager::DiscardObservation(bool discard_observation)
{
	discard_observation_ = discard_observation;
}

//! Returns the observations.
/*! This method is called after 'SetTime' set the time at which the
 observations are requested.
 \param[out] observation observation vector.
 */
void SimvascularObservationManager::GetObservation(
		SimvascularObservationManager::observation& observation) {
	throw ErrorUndefined(
			"void SimvascularObservationManager::GetObservation(const state& x,"
			"SimvascularObservationManager::observation& observation)");
}

////////////////
// INNOVATION //
////////////////

//! Returns an innovation.
/*! This method is called after 'SetTime' set the time at which the
 innovation is requested.
 \param[in] state state vector.
 \param[out] innovation innovation vector.
 */
template<class state>
void SimvascularObservationManager::GetInnovation(const state& x,
		SimvascularObservationManager::observation& innovation) {

	//
	// for now, load the full state from an existing set of results
	//

	int lower_bound;
	int upper_bound;

	lower_bound = (int)(Nskip_*floor((time_-initial_time_)/(double)Nskip_));
	upper_bound = lower_bound + Nskip_;

	//
	// only load data if the time interval has changed
	//

	if (lower_bound != current_lower_bound_) {
		cout << "loading data at time " << lower_bound << endl;
		loadrestart(lower_bound,soln_lower_,acc_lower_,disp_lower_);
		current_lower_bound_ = lower_bound;
	}

	if (upper_bound != current_upper_bound_) {
		if (upper_bound <= final_time_) {
			cout << "loading data at time " << upper_bound << endl;
			loadrestart(upper_bound,soln_upper_,acc_upper_,disp_upper_);
		}
		current_upper_bound_ = upper_bound;
	}

	//
    // now we need to actually compute the innovation
    // compute the interpolation factors
	//

    double t_alpha = (time_ - current_lower_bound_)/(current_upper_bound_-current_lower_bound_);
    t_alpha = 1 - t_alpha;

    //cout << "alpha parameter " << t_alpha << endl;

    int obsCounter = 0;

    innovation.Reallocate(Nobservation_);

    for (int kk = 0; kk < Nobservation_; kk++) {
    	innovation(kk) = -x(StateObsIndex_(kk)) +
    			t_alpha*dataarrays_lower_[DataArraysObsIndex_(kk)] + (1-t_alpha)*dataarrays_upper_[DataArraysObsIndex_(kk)];
    }

    //cout << "stiffness " << pow(2.0,x(47344)) << endl;

}

////////////
// ACCESS //
////////////

//! Indicates if some observations are available at a given time.
/*!
 \param[in] time a given time.
 */
bool SimvascularObservationManager::HasObservation(double time) {
	return time_ <= final_time_ && time_ >= initial_time_;
}

//! Indicates if some observations are available at current time.
bool SimvascularObservationManager::HasObservation() const {
	return time_ <= final_time_ && time_ >= initial_time_;
}

//! Returns the number of available observations.
/*!
 \return The total number of observation at current time.
 */
int SimvascularObservationManager::GetNobservation() const {
	return Nobservation_;
}

///////////////
// OPERATORS //
///////////////

//! Applies the observation operator to a given vector.
/*! This method is called after 'SetTime' set the time at which the
 operator is defined.
 \param[in] x a vector.
 \param[out] y the value of the operator applied to \a x. It is resized
 if needed.
 */
template<class state>
void SimvascularObservationManager::ApplyOperator(const state& x,
		observation& y) const {

	int obsCounter = 0;

	for (int jj = 0; jj < isize_nshg_; jj++) {
		for (int kk = 0; kk < 4; kk++) {
			if (linobs_soln_[kk*isize_nshg_+jj] > 0) {
				y(obsCounter) = linobs_soln_[kk*isize_nshg_+jj]*x(kk*isize_nshg_+jj);
				obsCounter++;
			}
		}
	}

	for (int jj = 0; jj < isize_nshg_; jj++) {
		for (int kk = 0; kk < 4; kk++) {
			if (linobs_acc_[kk*isize_nshg_+jj] > 0) {
				y(obsCounter) = linobs_acc_[kk*isize_nshg_+jj]*x(kk*isize_nshg_+jj + isize_nshg_*4);
				obsCounter++;
			}
		}
	}

	for (int jj = 0; jj < isize_nshg_; jj++) {
		for (int kk = 0; kk < 3; kk++) {
			if (linobs_disp_[kk*isize_nshg_+jj] > 0) {
				y(obsCounter) = linobs_disp_[kk*isize_nshg_+jj]*x(kk*isize_nshg_+jj + isize_nshg_*4 + isize_nshg_*4);
				obsCounter++;
			}
		}
	}

}

//! Return an observation error covariance.
/*!
 \param[in] i row index.
 \param[in] j column index.
 \return The element (\a i, \a j) of the observation error variance.
 */
double SimvascularObservationManager::GetErrorVariance(int i, int j) const {
	throw ErrorUndefined("double SimvascularObservationManager"
			"::GetErrorVariance(int i, int j) const");
}

//! Returns the observation error variance.
/*!
 \return The observation error covariance matrix.
 */
const SimvascularObservationManager::error_variance&
SimvascularObservationManager::GetErrorVariance() const {
	throw ErrorUndefined("const typename"
			"SimvascularObservationManager::error_variance& "
			"SimvascularObservationManager::GetErrorVariance() const");
}

//! Returns the inverse of the observation error covariance matrix.
/*!
 \return The inverse of the matrix of the observation error covariance.
 */
const SimvascularObservationManager::error_variance&
SimvascularObservationManager::GetErrorVarianceInverse() const {
	return error_variance_inverse_;
}

void SimvascularObservationManager::loadrestart(int timeindex, double* soln, double* acc, double* disp) {

	int irestart; /* file handle for restart */
	int iarray[10];
	int ithree = 3;
	int isize;
	char iformat[80];
	char filename[255];
	char restart_filename[255];

	strcpy(iformat, outpar.iotype);

	// assumed that the correct directory is already set
	strcpy(filename,data_directory_.c_str());
	sprintf(restart_filename, "restart.%d.1", timeindex);
	strcat(filename,"/");
	strcat(filename,restart_filename);

	cout << "trying to load " << filename << endl;

	openfile_(filename, "read", &irestart);

	// read the state arrays
	readheader_(&irestart, "solution?", (void*) iarray, &ithree, "double",
				iformat);
	isize = iarray[0]*iarray[1];
	readdatablock_(&irestart, "solution?", (void*) soln, &isize, "double",
			iformat);

	readheader_(&irestart, "time derivative of solution?", (void*) iarray, &ithree, "double",
					iformat);
	isize = iarray[0]*iarray[1];
	readdatablock_(&irestart, "time derivative of solution?", (void*) acc, &isize,
			"double", iformat);

	readheader_(&irestart, "displacement?", (void*) iarray, &ithree, "double",
			iformat);
	isize = iarray[0]*iarray[1];
	readdatablock_(&irestart, "displacement?", (void*) disp, &isize,
			"double", iformat);



}

//! Returns the name of the class.
/*!
 \return The name of the class.
 */
string SimvascularObservationManager::GetName() const {
	return "SimvascularObservationManager";
}

//! Receives and handles a message.
/*
 \param[in] message the received message.
 */
void SimvascularObservationManager::Message(string message) {
	// Put here any processing you need.
}

}

#define SIMVASCULAROBSERVATIONMANAGER_CXX
#endif
