#ifndef SIMVASCULAROBSERVATIONMANAGER_CXX

#include "SimvascularObservationManager.hxx"

namespace Verdandi {

/////////////////////////////////
// CONSTRUCTORS AND DESTRUCTOR //
/////////////////////////////////

//! Default constructor.
SimvascularObservationManager::SimvascularObservationManager()
:   Nobservation_(0),
    Nobservation_local_(0),
    Nobservation_nodal_(0),
    Nobservation_dist_(0),
    Nobservation_flow_(0),
    Nobservation_avgpressure_(0),
    Nskip_(0),
    initial_time_(0),
    final_time_(0),
    use_restarts_(0),
    execute_nodal_observations_(0),
    execute_distance_observations_(0),
    time_(0),
    data_period_(0),
    discard_observation_(0),
    error_variance_value_(1),
    error_variance_value_nodal_(1),
    error_variance_value_dist_(1),
    error_variance_value_avgpress_(1),
    error_variance_value_flow_(1),
    Nstate_model_(0),
    isize_solution_(0),
    isize_displacement_(0),
    isize_nshg_(0),
    isize_nshguniq_(0),
    rank_(0),
    numProcs_(1),
    current_lower_bound_(-1),
    current_upper_bound_(-1),
    gat(NULL)
{

}

//! Destructor.
SimvascularObservationManager::~SimvascularObservationManager() {
	// Operations to be performed when the object is destroyed.

	obs_out_part_.close();
	obs_out_single_.close();

	obs_in_part_.close();
	obs_out_single_.close();

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

	MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

	MPI_Comm_size(MPI_COMM_WORLD, &numProcs_);

	Nstate_model_ = model.GetNstate();
	Nobservation_local_ = 0;

    int obsCounter = 0;

	configuration.SetPrefix("observation.");
	configuration.Set("use_restarts",use_restarts_);
	configuration.Set("data_directory", data_directory_);
	configuration.Set("Nskip", "v > 0", Nskip_);
	configuration.Set("initial_time", "", 0., initial_time_);
	configuration.Set("final_time", "", numeric_limits<double>::max(),final_time_);
	configuration.Set("data_period","",final_time_,data_period_);

	configuration.Set("execute_nodal_observations_",execute_nodal_observations_);

	configuration.Set("execute_distance_observations_",execute_distance_observations_);

	configuration.Set("Nobservation_flow",Nobservation_flow_);

	configuration.Set("Nobservation_avgpressure",Nobservation_avgpressure_);

	if (Nobservation_flow_ > 0 || Nobservation_avgpressure_ > 0) {
		csobs_origins_.resize(Nobservation_flow_+Nobservation_avgpressure_);
		csobs_normals_.resize(Nobservation_flow_+Nobservation_avgpressure_);

		configuration.Set("csobs_origins",csobs_origins_);
		configuration.Set("csobs_normals",csobs_normals_);
		configuration.Set("csobs_radii",csobs_radii_);
	}

	configuration.Set("error.variance", "v > 0", error_variance_value_);

	configuration.Set("error.variance_nodal", "v > 0", error_variance_value_nodal_);
	configuration.Set("error.variance_dist", "v > 0", error_variance_value_dist_);
	cout << "error variance dist " << error_variance_value_dist_ << endl;
	configuration.Set("error.variance_avgpress", "v > 0", error_variance_value_avgpress_);
	cout << "error variance p avg " << error_variance_value_avgpress_ << endl;
	configuration.Set("error.variance_flow", "v > 0", error_variance_value_flow_);

	isize_solution_ = conpar.nshg * 4; // ignore last dof

	if (nomodule.ideformwall > 0)
		isize_displacement_ = conpar.nshg * NSD ;
	else
		isize_displacement_ = 0;

	isize_nshg_ = conpar.nshg;

	// Get pointer to the single instance of PhGlobalArrayTransfer
	gat = PhGlobalArrayTransfer::Instance();

	// Set up observation operator and increment number of local observations
    int actualnshg = conpar.nshguniq;

    // -------------------------------------------------------
    // 1. compact way of representing single node observations
    if (execute_nodal_observations_) {
    	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

    		int actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];
    		actualIdx--;

    		for (int kk = 0; kk < 4; kk++) {

    			int obsFuncVal = (gat->global_ilinobsfunc_sol_ptr)[kk * conpar.nshg + actualIdx];

    			if (obsFuncVal > 0) {
    				StateObsIndex_.PushBack(kk+4*unitIdx);
    				//if (use_restarts_)
    				//	DataArraysObsIndex_.PushBack(kk*isize_nshg_+actualIdx);
    				//else
    				DataArraysObsIndex_.PushBack(obsCounter);
    				obsCounter++;
    			}
    		}
    	}

    	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

    		int actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];
    		actualIdx--;

    		for (int kk = 0; kk < 4; kk++) {

    			int obsFuncVal = (gat->global_ilinobsfunc_acc_ptr)[kk * conpar.nshg + actualIdx];

    			if (obsFuncVal > 0) {
    				StateObsIndex_.PushBack(kk+4*unitIdx + actualnshg*4);
    				//if (use_restarts_)
    				//	DataArraysObsIndex_.PushBack(kk*isize_nshg_+actualIdx + isize_solution_);
    				//else
    				DataArraysObsIndex_.PushBack(obsCounter);
    				obsCounter++;
    			}
    		}
    	}

    	if (nomodule.ideformwall > 0) {

    		for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

    			int actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];
    			actualIdx--;

    			for (int kk = 0; kk < 3; kk++) {

    				int obsFuncVal = (gat->global_ilinobsfunc_disp_ptr)[kk * conpar.nshg + actualIdx];

    				if (obsFuncVal > 0) {
    					StateObsIndex_.PushBack(kk+3*unitIdx + actualnshg*4 + actualnshg*4);
    					//if (use_restarts_)
    					//	DataArraysObsIndex_.PushBack(kk*isize_nshg_+actualIdx + isize_solution_*2);
    					//else
    					DataArraysObsIndex_.PushBack(obsCounter);
    					obsCounter++;
    				}
    			}
    		}
    	}
    }

    Nobservation_nodal_ = obsCounter;
    Nobservation_local_ += Nobservation_nodal_;

    // -------------------------------
    // 2. set up distance observation

    // the measured distance will be read from the global arrays and will be used
    // directly in getinnovation
    // here, we simply increment the local observation count

    //obsCounter = 0;
    Nobservation_dist_ = 0;

    if (execute_distance_observations_ && nomodule.imeasdist > 0) {

    	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

    		int actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];
    		actualIdx--;

    		int obsFuncVal = (gat->global_obsfunc_dist_ptr)[actualIdx];

    		if (obsFuncVal > 0) {
    			// this is the position in the state vector of the distance vector
    			StateObsIndex_.PushBack(unitIdx + actualnshg*3 + actualnshg*4 + actualnshg*4);
    			DataArraysObsIndex_.PushBack(obsCounter);
    			obsCounter++;
    			Nobservation_dist_++;
    		}
    	}

    }

    Nobservation_local_ += Nobservation_dist_;

    // -----------------------------------------------------
    // 3. set up cross-sectional flow and pressure observation with VTK

    if (Nobservation_flow_+Nobservation_avgpressure_ > 0) {

    	geom_points_ = vtkSmartPointer<vtkPoints>::New();
    	geom_ids_ = vtkSmartPointer<vtkIdList>::New();
    	geom_UGrid_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
    	vtkSmartPointer<vtkDoubleArray> geom_vel_array =
    			vtkSmartPointer<vtkDoubleArray>::New();
    	vtkSmartPointer<vtkDoubleArray> geom_pres_array =
    			vtkSmartPointer<vtkDoubleArray>::New();

    	geom_vel_array->SetNumberOfComponents(3);

    	for (int unitIdx = 0; unitIdx < conpar.nshg; unitIdx++) {

    		double coordVal[3];
    		double solVal[4];

    		for (int kk = 0; kk < 3; kk++)
    			coordVal[kk] = (gat->global_coord_ptr)[kk * conpar.nshg + unitIdx];

    		for (int kk = 0; kk < 4; kk++)
    			solVal[kk] = (gat->global_yold_ptr)[kk * conpar.nshg + unitIdx];

    		geom_points_->InsertPoint(unitIdx,coordVal[0],coordVal[1],coordVal[2]);
    		geom_vel_array->InsertNextTuple3(solVal[0],solVal[1],solVal[2]);
    		geom_pres_array->InsertNextTuple1(solVal[3]);
    	}

    	geom_UGrid_->SetPoints(geom_points_);

    	for (int kk = 0; kk < gat->global_mien.size(); kk++)
    		for (int jj = 0; jj < gat->global_npro[kk]; jj++) {

    			for (int ii = 0; ii < 4; ii++) {
    				int nodeIndex = gat->global_mien[kk][ii * gat->global_npro[kk] + jj];
    				geom_ids_->InsertNextId(--nodeIndex);
    			}

    			geom_UGrid_->InsertNextCell(VTK_TETRA,geom_ids_);

    			geom_ids_->Reset();
    		}

    	geom_UGrid_->GetPointData()->AddArray(geom_vel_array);
    	geom_UGrid_->GetPointData()->GetArray(0)->SetName("velocity");
    	geom_UGrid_->GetPointData()->AddArray(geom_pres_array);
    	geom_UGrid_->GetPointData()->GetArray(1)->SetName("pressure");
    	geom_UGrid_->Update();

    	// here we specify the number of cross-sections we want to observe flow on
    	// the single slice plane may cut through several parts of the
    	// mesh and only one part (i.e. the local geometry) is the desired cut.
    	// the desired cut is not necssarily connected and may not be on the current processor

    	for (int kk = 0; kk < Nobservation_flow_+Nobservation_avgpressure_; kk++) {

    		// set the location of the cutting plane
    		geom_plane_ = vtkSmartPointer<vtkPlane>::New();
    		geom_plane_->SetOrigin(csobs_origins_[kk](0),csobs_origins_[kk](1),csobs_origins_[kk](2));
    		geom_plane_->SetNormal(csobs_normals_[kk](0),csobs_normals_[kk](1),csobs_normals_[kk](2));
    		geom_planes_.push_back(geom_plane_);

    		geom_cutter_ = vtkSmartPointer<vtkCutter>::New();
    		geom_cutter_->SetCutFunction(geom_plane_);
    		geom_cutter_->SetInput(geom_UGrid_);
    		geom_cutter_->Update();

    		geom_connectivity_ = vtkSmartPointer<vtkConnectivityFilter>::New();
    		geom_connectivity_->SetInputConnection(geom_cutter_->GetOutputPort());
    		geom_connectivity_->SetExtractionModeToClosestPointRegion();
    		geom_connectivity_->SetClosestPoint(csobs_origins_[kk](0),csobs_origins_[kk](1),csobs_origins_[kk](2));
    		geom_connectivity_->Update();

    		vector <double> distfromorigin;
    		double testpoint[3];

    		// measure the distance from the origin to each cell of the cut
    		// if the distance is greater than a user specified threshold,
    		// ignore that cell in the calculation
    		for (int jj = 0; jj < 3; jj++)
    			testpoint[jj] = csobs_origins_[kk](jj);

    		for (int jj = 0; jj < geom_cutter_->GetOutput()->GetNumberOfCells(); jj++) {

    			double coord[3][3],closestpoint[3];
    			int locinfo;
    			double *tempcoord;
    			vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

    			geom_cutter_->GetOutput()->GetCellPoints(jj,ptIds);

    			for (int ii = 0; ii < 3; ii++) {

    				tempcoord = geom_cutter_->GetOutput()->GetPoint(ptIds->GetId(ii));
    				for (int nn = 0; nn < 3; nn++)
    					coord[ii][nn] = tempcoord[nn];

    			}

    			cptrip( &testpoint[0], &coord[0][0], &coord[1][0], &coord[2][0], &closestpoint[0], locinfo );

    			double closestdistsqr = (closestpoint[0]-testpoint[0])*(closestpoint[0]-testpoint[0])+
    					(closestpoint[1]-testpoint[1])*(closestpoint[1]-testpoint[1])+
    					(closestpoint[2]-testpoint[2])*(closestpoint[2]-testpoint[2]);

    			distfromorigin.push_back(sqrt(closestdistsqr));

    		}

    		distances_fromorigin_.push_back(distfromorigin);

    		geom_cutters_.push_back(geom_cutter_);
    		geom_connec_filters_.push_back(geom_connectivity_);

    		//        vtkSmartPointer<vtkPolyDataWriter> writer = vtkPolyDataWriter::New();
    		//        ostream *vtkout;
    		//
    		//        if (rank_ == 0)
    		//        	writer->SetFileName("fromcutter_1.vtk");
    		//        else
    		//        	writer->SetFileName("fromcutter_2.vtk");
    		//        writer->SetInput(geom_cutter_->GetOutput());
    		//        vtkout = writer->OpenVTKFile();
    		//        writer->Write();
    		//        writer->CloseVTKFile(vtkout);

    		// number of local observations
    		if (rank_ == numProcs_ - 1) {
    			DataArraysObsIndex_.PushBack(obsCounter++);

    			// increment local observation count
    			Nobservation_local_++;
    		}

    	}

    }

    // compute the global number of observation
    MPI_Allreduce(&Nobservation_local_, &Nobservation_, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	// we don't need to populate data arrays for distance observations
	// distance data is handled directly in getinnovation
	dataarrays_lower_.Reallocate(Nobservation_local_);
	dataarrays_upper_.Reallocate(Nobservation_local_);
	dataarrays_lower_.Zero();
	dataarrays_upper_.Zero();

	// observation error variance matrix
#ifdef VERDANDI_ROUKF_PARALLEL_INNOVATION

	// here we assume that is diagonal
	error_variance_inverse_diag_.Reallocate(Nobservation_, Nobservation_local_);
	int local_start, local_end;
	error_variance_inverse_diag_.GetProcessorRange(local_start, local_end);

	// TODO
//	for (int kk = 0; kk < Nobservation_local_; kk++) {
//		error_variance_inverse_diag_.SetBuffer( local_start+kk, double(1) / error_variance_value_ );
//	}
	int ncounter = 0;

	for (int kk = 0; kk < Nobservation_nodal_; kk++) {
		error_variance_inverse_diag_.SetBuffer( local_start+ncounter, double(1) / error_variance_value_nodal_ );
		ncounter++;
	}
	for (int kk = 0; kk < Nobservation_dist_; kk++) {
		error_variance_inverse_diag_.SetBuffer( local_start+ncounter, double(1) / error_variance_value_dist_ );
		ncounter++;
	}
	if (rank_ == numProcs_ - 1) {
		for (int kk = 0; kk < Nobservation_flow_; kk++) {
			error_variance_inverse_diag_.SetBuffer( local_start+ncounter, double(1) / error_variance_value_flow_ );
			ncounter++;
		}
		for (int kk = 0; kk < Nobservation_avgpressure_; kk++) {
			error_variance_inverse_diag_.SetBuffer( local_start+ncounter, double(1) / error_variance_value_avgpress_ );
			ncounter++;
		}
	}
	error_variance_inverse_diag_.Flush();

#else

#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
	build_diagonal_sparse_matrix(Nobservation_, error_variance_value_,
			error_variance_);
	build_diagonal_sparse_matrix(Nobservation_,
			double(double(1) / error_variance_value_),
			error_variance_inverse_);
#else
	error_variance_.Reallocate(Nobservation_, Nobservation_);
	error_variance_.SetIdentity();
	Mlt(error_variance_value_, error_variance_);
	error_variance_inverse_.Reallocate(Nobservation_, Nobservation_);
	error_variance_inverse_.SetIdentity();
	Mlt(double(double(1)/ error_variance_value_), error_variance_inverse_);
#endif

#endif

    // some console output
	if (rank_ == 0)
		cout << "Simvascular Observation Manager initiated" << endl;

	cout << "At rank " << rank_ << " we have " << Nobservation_local_ << " observations ";
	cout << Nobservation_nodal_ << " nodal ";
	cout << Nobservation_dist_ << " dist " << endl << endl;

	if (rank_ == 0) {
		cout << "global flow observations " << Nobservation_flow_ << endl;
		cout << "global avgpress observations " << Nobservation_avgpressure_ << endl;
	}

	// files for output (only used for generating synthetic data)
	stringstream s_temp;
	s_temp << rank_;
	obsfilename_part_ = "saved_observations." + s_temp.str() + ".dat";
	obsfilename_single_ = "saved_observations_single.dat";

	obs_out_part_.open (obsfilename_part_.c_str(), ios::out | ios::app );
	obs_out_single_.open (obsfilename_single_.c_str(), ios::out | ios::app );

	obs_out_part_ << std::setprecision( std::numeric_limits<double>::digits10+1);
	obs_out_single_ << std::setprecision( std::numeric_limits<double>::digits10+1);

    obsfilename_part_ = data_directory_ + "/" + obsfilename_part_;
    obsfilename_single_ = data_directory_ + "/" + obsfilename_single_;

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
void SimvascularObservationManager::GetInnovation(const state& x, observation& innovation) {
	throw ErrorUndefined(
			"void SimvascularObservationManager::"
			"GetInnovation(const state& x, observation& innovation)");
}

//! Returns an innovation.
/*! This method is called after 'SetTime' set the time at which the
 innovation is requested.
 \param[in] state state vector.
 \param[out] innovation innovation vector.
 */
#if defined(VERDANDI_ROUKF_PARALLEL_INNOVATION)
template<class state>
void SimvascularObservationManager::GetInnovation(const state& x,
		state& innovation_p_orig,
		state& innovation_p_fe) {
#else
template<class state>
void SimvascularObservationManager::GetInnovation(const state& x,
		SimvascularObservationManager::observation& innovation_orig,
		SimvascularObservationManager::observation& innovation_fe) {
#endif

	// for now, load the full state from an existing set of results
	int lower_bound;
	int upper_bound;
	int period_time;

	period_time = (int)time_ % (int)data_period_;

	lower_bound = (int)(Nskip_*floor(period_time/(double)Nskip_));
	upper_bound = lower_bound + Nskip_;

	// only load data if the time interval has changed
	if ((int)time_ < final_time_) {
		if (lower_bound != current_lower_bound_) {
			if (rank_ == 0)
				cout << "loading data at time " << lower_bound << endl;

			//if (use_restarts_)
			//	loadrestart(lower_bound, soln_lower_, acc_lower_, disp_lower_);
			//else
			LoadObservationSingleLocal(lower_bound, dataarrays_lower_);

			current_lower_bound_ = lower_bound;
		}

		if (upper_bound != current_upper_bound_) {

			if (rank_ == 0)
				cout << "loading data at time " << upper_bound << endl;

			//if (use_restarts_)
			//	loadrestart(upper_bound, soln_upper_, acc_upper_, disp_upper_);
			//else
			LoadObservationSingleLocal(upper_bound, dataarrays_upper_);

			current_upper_bound_ = upper_bound;
		}
	}

    // compute the interpolation factors
    double t_alpha = ((int)time_ - current_lower_bound_)/(current_upper_bound_-current_lower_bound_);
    t_alpha = 1 - t_alpha;

    observation zHx1,zHx2;

//    // debugging output
//    std::ostringstream ostr; //output string stream
//    std::string filename = "zHx";
//    ostr << rank_;
//    filename = filename+ostr.str();
//    ofstream outfile;
//	outfile.open(filename.c_str());

    this->ApplyOperatorLocal(x,zHx1,zHx2);

    for (int kk = 0; kk < Nobservation_local_; kk++) {

    	zHx1(kk) = -zHx1(kk) +
    			t_alpha*dataarrays_lower_(DataArraysObsIndex_(kk)) +
    			(1-t_alpha)*dataarrays_upper_(DataArraysObsIndex_(kk)) ;

    	zHx2(kk) = -zHx2(kk) +
    			t_alpha*dataarrays_lower_(DataArraysObsIndex_(kk)) +
    			(1-t_alpha)*dataarrays_upper_(DataArraysObsIndex_(kk)) ;

//    	cout << "zhx1 " << zHx1(kk) << " zhx2 " << zHx2(kk) << endl;
//    	cout << " l " << dataarrays_lower_(DataArraysObsIndex_(kk));
//    	cout << " u " << dataarrays_upper_(DataArraysObsIndex_(kk)) << endl;

    }

    //outfile.close();
    //

#if defined(VERDANDI_ROUKF_PARALLEL_INNOVATION)
    // distributed innovation vector
    innovation_p_orig.Reallocate(Nobservation_, Nobservation_local_);
    innovation_p_fe.Reallocate(Nobservation_, Nobservation_local_);

    int local_start, local_end;
    innovation_p_orig.GetProcessorRange(local_start, local_end);

    for (int kk = 0; kk < Nobservation_local_; kk++) {
    	innovation_p_orig.SetBuffer( local_start+kk ,zHx1(kk) );

    	innovation_p_fe.SetBuffer( local_start+kk, error_variance_inverse_diag_(local_start+kk) * zHx2(kk) );
    }
    innovation_p_orig.Flush();
    innovation_p_fe.Flush();

#else
    // the innovation is required to be on all processors
    innovation_orig.Reallocate(Nobservation_);
    innovation_orig.Zero();

    innovation_fe.Reallocate(Nobservation_);
    innovation_fe.Zero();

    int *obs_recvcount = new int[numProcs_];
    int *obs_displ = new int[numProcs_];

    MPI_Allgather(&Nobservation_local_, 1, MPI_INT,
    		obs_recvcount, 1, MPI_INT,
    		MPI_COMM_WORLD);

    obs_displ[0] = 0;
    for (int kk = 1; kk < numProcs_; kk++) {
    	obs_displ[kk] = obs_displ[kk-1]+obs_recvcount[kk-1];
    }

    //	for (int kk = 0; kk < numProcs_ ; kk++) {
    //		cout << "receive count " << obs_recvcount[kk] << " ";
    //	}
    //	cout << endl;

    MPI_Allgatherv(zHx1.GetData(), Nobservation_local_, MPI_DOUBLE,
                   innovation_orig.GetData(), obs_recvcount, obs_displ, MPI_DOUBLE,
                   MPI_COMM_WORLD);

    MPI_Allgatherv(zHx2.GetData(), Nobservation_local_, MPI_DOUBLE,
    		       innovation_fe.GetData(), obs_recvcount, obs_displ, MPI_DOUBLE,
    		       MPI_COMM_WORLD);

    delete [] obs_recvcount;
    delete [] obs_displ;

#endif
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

//! Returns the size of the local (on-processor) state vector.
/*!
      \return The size of the local state vector.
 */
int SimvascularObservationManager::GetLocalNobservation() const {
	return Nobservation_local_; // this is used in reallocate routine for petsc matrices in ROUKF
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
void SimvascularObservationManager::ApplyOperator(const state& x, observation& y) const {

	// not implemented as of now since we will always use the getInnovation

	throw ErrorUndefined(
				"void SimvascularObservationManager::ApplyOperator(const state& x, observation& y)");

}

//! Applies the observation operator to a given vector.
/*! This method is called after 'SetTime' set the time at which the
 operator is defined.
 \param[in] x a vector.
 \param[out] y the value of the operator applied to \a x. It is resized
 if needed.
 */
template<class state>
void SimvascularObservationManager::ApplyOperatorLocal(const state& x, observation& Hx1, observation& Hx2) {

	observation Hx_cs;

	int icounter = 0;
	int ncounter = 0;

	int state_start, state_end;
	x.GetProcessorRange(state_start, state_end);

	Hx1.Reallocate(Nobservation_local_);
	Hx2.Reallocate(Nobservation_local_);

	Hx1.Zero();
	Hx2.Zero();

	// simple nodal observation
	for (int kk = 0; kk < Nobservation_nodal_; kk++) {

		Hx1(kk) = x(state_start+StateObsIndex_(icounter));
		Hx2(kk) = x(state_start+StateObsIndex_(icounter));

		icounter++;
	}

	// distance observations
	ncounter += icounter;
	icounter = 0;
	int actualnshg = conpar.nshguniq;

	for (int kk = 0; kk < Nobservation_dist_; kk++) {

		// grab the distance vector
		Hx1(ncounter+kk) = x(state_start+StateObsIndex_(ncounter+icounter));

		//skip ahead to grab the finite element distance vector
		Hx2(ncounter+kk) = x(state_start+StateObsIndex_(ncounter+icounter)+actualnshg);

		icounter++;
	}

	// flow observation
	// the values are located on the last processor
	ncounter += icounter;
	icounter = 0;

	if (Nobservation_flow_ + Nobservation_avgpressure_ > 0)
		ApplyOperatorFlow(x,Hx_cs);

	if (rank_ == numProcs_ - 1) {
		for (int kk = 0; kk < Nobservation_flow_ + Nobservation_avgpressure_; kk++) {

			Hx1(ncounter+kk) = Hx_cs(kk);
			Hx2(ncounter+kk) = Hx_cs(kk);

			//cout << Hx_start+ncounter+kk << " " << Hx_flow(kk) << endl;

			icounter++;
		}
	}

}

template<class state>
void SimvascularObservationManager::ApplyOperatorFlow(const state& x, observation& Hx) {

	int state_start, state_end, icounter;

    x.GetProcessorRange(state_start, state_end);

    icounter = state_start;

	vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

	// Here we look through the state vector and grab the velocity values
    // to assign to a temporary array that will allow communication of those values to the
	// duplicated nodes on the interprocessor boundaries

	// the use of this temporary array could possibly
	// be avoided if model stateupdated is called
	// before the getinnovation function is called
	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

		int actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond

			double val = x(icounter++);

			(gat->global_temporary_array_ptr)[varIdx * conpar.nshg + actualIdx-1] = val;

		}

	}

	// the values coming from the state vector need to communicated at the interprocessor boundaries
	estim_helpers_temp_comm();

	// now we reassign the velocities values to update the cross-sectional avg flow

	for (int unitIdx = 0; unitIdx < conpar.nshg; unitIdx++) {

		double solVal[4];

		for (int ii = 0; ii < 4; ii++)
			solVal[ii] = (gat->global_temporary_array_ptr)[ii * conpar.nshg + unitIdx];

		geom_UGrid_->GetPointData()->GetArray(0)->SetTuple3(unitIdx,solVal[0],solVal[1],solVal[2]);
		geom_UGrid_->GetPointData()->GetArray(1)->SetTuple1(unitIdx,solVal[3]);

	}

	if (rank_ == numProcs_ - 1)
		Hx.Reallocate(Nobservation_flow_+Nobservation_avgpressure_);

	for (int kk = 0; kk < Nobservation_flow_ + Nobservation_avgpressure_; kk++) {

		double area_local, area, avgFlow_local, avgFlow, avgPres_local, avgPres;

		// the next two calls update the values on cross-sectional cut
		geom_cutters_[kk]->Modified();
		geom_cutters_[kk]->Update();
		//geom_connec_filters_[kk]->Update();

		//	for (int kk = 0; kk < geom_connectivity_->GetOutput()->GetNumberOfPoints(); kk++) {
		//		coord1 = geom_connectivity_->GetOutput()->GetPoint(kk);
		//		cout << kk << " " << coord1[0] << " " << coord1[1] << " " << coord1[2] << endl;
		//	}

		// loop through cells of the cut to compute
		// the cross-sectional avg flow

		avgFlow_local = 0;
		area_local = 0;
		avgPres_local = 0;
		for (int jj = 0; jj < geom_cutters_[kk]->GetOutput()->GetNumberOfCells(); jj++) {

			if (distances_fromorigin_[kk][jj] <= csobs_radii_[kk]) {

				double coord[3][3];
				double vel[3][3],press123[3];
				double *tempcoord,*tempvel;
				double A[3], B[3], C[3];
				double triArea;
				double tempL;

				geom_cutters_[kk]->GetOutput()->GetCellPoints(jj,ptIds);


				for (int ii = 0; ii < 3; ii++) {

					tempvel = geom_cutters_[kk]->GetOutput()->GetPointData()->GetArray(0)->GetTuple3(ptIds->GetId(ii));

					tempcoord = geom_cutters_[kk]->GetOutput()->GetPoint(ptIds->GetId(ii));

					for (int nn = 0; nn < 3; nn++)
						vel[ii][nn] = tempvel[nn];

					for (int nn = 0; nn < 3; nn++)
						coord[ii][nn] = tempcoord[nn];

				}

				for (int ii = 0; ii < 3; ii++)
					press123[ii] = geom_cutters_[kk]->GetOutput()->GetPointData()->GetArray(1)->GetTuple1(ptIds->GetId(ii));

				for (int ii = 0; ii < 3; ii++) {
					A[ii] = coord[1][ii]-coord[0][ii];
					B[ii] = coord[2][ii]-coord[0][ii];
				}

				C[0] = (A[1]*B[2])-(B[1]*A[2]);
				C[1] = -(A[0]*B[2])+(B[0]*A[2]);
				C[2] = (A[0]*B[1])-(A[1]*B[0]);

				tempL = sqrt(C[0]*C[0]+C[1]*C[1]+C[2]*C[2]);

				triArea = 0.5*tempL;

				for (int ii = 0; ii < 3; ii++)
					C[ii] = C[ii] / tempL;

				area_local += triArea;

				avgFlow_local += (double(1.0)/3.0)*(vel[0][0]*C[0]+vel[0][1]*C[1]+vel[0][2]*C[2]+
						vel[1][0]*C[0]+vel[1][1]*C[1]+vel[1][2]*C[2]+
						vel[2][0]*C[0]+vel[2][1]*C[1]+vel[2][2]*C[2])*triArea;

				avgPres_local += (double(1.0)/3.0)*(press123[0]+press123[1]+press123[2])*triArea;
			}

		}

		// Compute the total flow and average pressure by summing across processes
		MPI_Reduce(&avgFlow_local, &avgFlow, 1, MPI_DOUBLE, MPI_SUM, numProcs_ - 1,MPI_COMM_WORLD);
		MPI_Reduce(&avgPres_local, &avgPres, 1, MPI_DOUBLE, MPI_SUM, numProcs_ - 1,MPI_COMM_WORLD);
		MPI_Reduce(&area_local, &area, 1, MPI_DOUBLE, MPI_SUM, numProcs_ - 1,MPI_COMM_WORLD);

		if (rank_ == numProcs_ -1) {
			(kk < Nobservation_flow_) ?
				Hx(kk) = avgFlow : Hx(kk) = avgPres / area;
		}
	}



//	if (rank_ == numProcs_ -1)
//		this->flow_out_ << endl;
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

//void SimvascularObservationManager::loadrestart(int timeindex, double* soln, double* acc, double* disp) {
//
//	int irestart; /* file handle for restart */
//	int iarray[10];
//	int ithree = 3;
//	int isize;
//	char iformat[80];
//	char filename[255];
//	char restart_filename[255];
//
//	strcpy(iformat, outpar.iotype);
//
//	// assumed that the correct directory is already set
//	strcpy(filename,data_directory_.c_str());
//	sprintf(restart_filename, "restart.%d.%d", timeindex, rank_+1);
//	strcat(filename,"/");
//	strcat(filename,restart_filename);
//
//	//cout << "trying to load " << filename << endl;
//
//	openfile_(filename, "read", &irestart);
//
//	// read the state arrays
//	readheader_(&irestart, "solution?", (void*) iarray, &ithree, "double",
//				iformat);
//	isize = iarray[0]*iarray[1];
//	readdatablock_(&irestart, "solution?", (void*) soln, &isize, "double",
//			iformat);
//
//	readheader_(&irestart, "time derivative of solution?", (void*) iarray, &ithree, "double",
//			iformat);
//	isize = iarray[0]*iarray[1];
//	readdatablock_(&irestart, "time derivative of solution?", (void*) acc, &isize,
//			"double", iformat);
//
//	if (nomodule.ideformwall > 0) {
//		readheader_(&irestart, "displacement?", (void*) iarray, &ithree, "double",
//				iformat);
//		isize = iarray[0]*iarray[1];
//		readdatablock_(&irestart, "displacement?", (void*) disp, &isize,
//				"double", iformat);
//	}
//}

void SimvascularObservationManager::LoadObservationSingleLocal(int timeindex, Vector<double>& dataarray) {

	// in simvascular we assume that the time index is always an integral value
	// we need to sequentially read the file until
	int linetoread = timeindex / Nskip_;
	int icounter = 0;
	int ncounter = 0;
	string line;

    obs_in_part_.open(obsfilename_part_.c_str());

    // read in simple nodal observation data
    if (execute_nodal_observations_) {
    	if (obs_in_part_.is_open()) {

    		// skip to the desired line
    		while ( icounter < linetoread) {
    			getline (obs_in_part_,line);
    			icounter++;
    		}

    		// read in the values
    		for (int kk = 0; kk < Nobservation_nodal_; kk++) {
    			obs_in_part_ >> dataarray(ncounter);
    			ncounter++;
    		}

    		obs_in_part_.close();
    	}
    	else cout << "Unable to open file: " << obsfilename_part_ << endl;
    }

    if (execute_distance_observations_) {
    	if (obs_in_part_.is_open()) {

    		// simply increment counter since the
    		// distances are not loaded from this file
    		for (int kk = 0; kk < Nobservation_dist_; kk++)
    			ncounter++;
    	}
    }

    // read in flow observation data
    if (rank_ == numProcs_ - 1) {

    	obs_in_single_.open(obsfilename_single_.c_str());
    	icounter = 0;

    	if (obs_in_single_.is_open()) {

    		// skip to the desired line
    		while ( icounter < linetoread) {
    			getline (obs_in_single_,line);
    			icounter++;
    		}

    		// read in the values
    		for (int kk = 0; kk < Nobservation_flow_ + Nobservation_avgpressure_; kk++) {
    			obs_in_single_ >> dataarray(ncounter+kk);
    		}

    		obs_in_single_.close();
    	}

    	else cout << "Unable to open file: " << obsfilename_single_ << endl;
    	obs_in_single_.close();
    }

    obs_in_part_.close();

}

//! Saves observations in file
/*! This function isn't actually used in any data assimilator routine

 */
template<class state>
void SimvascularObservationManager::SaveObservationSingleLocal(const state& x) {

	// we assumed that the time has been set from the model time
	if ((int)time_ % Nskip_ == 0) {

		observation Hx1,Hx2;
		int ncounter = 0, ncounter2 = 0;

		cout << "SAVING OBSERVATIONS";

		this->ApplyOperatorLocal(x,Hx1,Hx2);

		//for (int kk = Hx_start; kk < Hx_end; kk++) {
		for (int kk = 0; kk < Nobservation_nodal_; kk++) {
			obs_out_part_ << Hx1(ncounter) << " ";
			ncounter++;
		}

		ncounter2 = ncounter;
		for (int kk = 0; kk < Nobservation_dist_; kk++) {
			obs_out_part_ << Hx1(ncounter) << " ";
			ncounter++;
		}

		for (int kk = 0; kk < Nobservation_dist_; kk++) {
			obs_out_part_ << Hx2(ncounter2) << " ";
			ncounter2++;
		}

		obs_out_part_ << endl;

		if (rank_ == numProcs_ - 1) {
			for (int kk = 0; kk < Nobservation_flow_+Nobservation_avgpressure_; kk++) {
				obs_out_single_ << Hx1(ncounter) << " ";
				ncounter++;
			}

			obs_out_single_ << endl;

		}

		cout << " [done] " << endl;

	}

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
