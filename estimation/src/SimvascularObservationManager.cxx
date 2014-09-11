#ifndef SIMVASCULAROBSERVATIONMANAGER_CXX
#define SIMVASCULAROBSERVATIONMANAGER_CXX

#include "SimvascularObservationManager.hxx"

namespace Verdandi {

/* */
/* */
SimvascularObservationType::SimvascularObservationType()
:Nobservation_local_(0),
 load_data_from_datfile_(0),
 isDistributed_(0)
{
	iNewComm_C_ = MPI_Comm_f2c(newcom.iNewComm);
	MPI_Comm_rank(iNewComm_C_, &rank_);
	MPI_Comm_size(iNewComm_C_, &numProcs_);
}


SimvascularObservationType::~SimvascularObservationType() {

}

void SimvascularObservationType::LoadData(int linetoread, Vector<double>& dataarray, int obs_start_index) {

	if (load_data_from_datfile_) {
		int icounter = 0;
		std::string line;

		data_file_in_.open(data_file_name_.c_str());

		if (data_file_in_.is_open()) {

			// skip to the desired line
			while ( icounter < linetoread) {
				std::getline (data_file_in_,line);
				icounter++;
			}

			// read in the values
			for (int kk = 0; kk < Nobservation_local_; kk++) {
				data_file_in_ >> dataarray(kk+obs_start_index);
			}

			data_file_in_.close();
		}
		else {
			throw Error("SimvascularObservationType::LoadData: Unable to open data file "+data_file_name_);
		}
	}

}

void SimvascularObservationType::SaveObservations(const state& x) {

	observation Hx1,Hx2;

	Hx1.Reallocate(Nobservation_local_);
	Hx2.Reallocate(Nobservation_local_);

	Hx1.Zero();
	Hx2.Zero();

	ApplyOperator(x,Hx1,Hx2,0);

	if (isDistributed_ || rank_ == numProcs_ - 1) {
		data_file_out_.open(saved_obs_file_name_.c_str(), ios::out | ios::app );
		data_file_out_ << std::setprecision( std::numeric_limits<double>::digits10+1);

		for (int kk = 0; kk < Nobservation_local_; ++kk)
			data_file_out_ << Hx1(kk) << " ";

		data_file_out_ << std::endl;
		data_file_out_.close();
	}

}

int SimvascularObservationType::getNobservation_local() const {
	return Nobservation_local_;
}

double SimvascularObservationType::getErrorVarianceValue(int ind) const {
	return error_variance_value_[ind];
}

std::string SimvascularObservationType::getName() const {
	return name_;
}

/* */
/* */
SimvascularNodalSolutionObservation::SimvascularNodalSolutionObservation() {
}

SimvascularNodalSolutionObservation::~SimvascularNodalSolutionObservation() {
}

void SimvascularNodalSolutionObservation::Initialize(std::string name, const SimvascularAugStatePart &aug_state_part, VerdandiOps &configuration) {

	int icounter = 0, ncounter = 0;

	name_ = name;
	std::cout << "Initializing " << name << std::endl;

	configuration.Set("data_directory", data_directory_);

	isDistributed_ = true;

	double error_variance_value_nodal;
	configuration.Set("error.variance_nodal", "v > 0", error_variance_value_nodal);

	// Get pointer to the single instance of SimvascularGlobalArrayTransfer
    gat = SimvascularGlobalArrayTransfer::Get();

	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {
		int actualIdx = (gat->pointerMapInt_["local index of unique nodes"])[unitIdx];
		actualIdx--;
		for (int kk = 0; kk < 4; kk++) {
			int obsFuncVal = (gat->pointerMapInt_["observation function solution"])[kk * conpar.nshg + actualIdx];
			if (obsFuncVal > 0) { // if we want to observe on this node
				state_obs_index_.push_back(aug_state_part.getDuplicatedStateIndex(ncounter));
				error_variance_value_.push_back(error_variance_value_nodal);
				icounter++;
			}
			ncounter++;
		}
	}

	Nobservation_local_ = icounter;

//	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {
//		int actualIdx = (gat->pointerMapInt_["local index of unique nodes"])[unitIdx];
//		actualIdx--;
//		for (int kk = 0; kk < 4; kk++) {
//			int obsFuncVal = (gat->pointerMapInt_["observation function time derivative of solution"])[kk * conpar.nshg + actualIdx];
//			if (obsFuncVal > 0) {
//				//StateObsIndex_.PushBack(kk+4*unitIdx + actualnshg*4);
//
//				//DataArraysObsIndex_.PushBack(obsCounter++);
//			}
//		}
//	}
//
	load_data_from_datfile_ = 1;
	data_file_name_ = "saved_observations_local_nodalsol";
	std::stringstream s_temp;
	s_temp << rank_;
	data_file_name_ = data_directory_ + "/" + data_file_name_ + s_temp.str() + ".dat";
	saved_obs_file_name_ = "saved_observations_local_nodalsol" + s_temp.str() + ".dat";

	// some console output
	std::cout << "[" << name_ << "] at rank (" << rank_ << ")" << " # of obs.: " << Nobservation_local_ << std::endl;

}

void SimvascularNodalSolutionObservation::ApplyOperator(const state& x, observation& Hx1, observation& Hx2, int obs_start_index) const {

	int icounter = 0;

	for (int kk = 0; kk < Nobservation_local_; kk++) {

		Hx1(obs_start_index+kk) = x(state_obs_index_[icounter]);
		Hx2(obs_start_index+kk) = Hx1(obs_start_index+kk);

		icounter++;
	}

}


/* */
/* */
SimvascularNodalDisplacementObservation::SimvascularNodalDisplacementObservation() {

}

SimvascularNodalDisplacementObservation::~SimvascularNodalDisplacementObservation() {

}

void SimvascularNodalDisplacementObservation::Initialize(std::string name, const SimvascularAugStatePart &aug_state_part, VerdandiOps &configuration) {

	int icounter = 0, ncounter = 0;

	name_ = name;
	std::cout << "Initializing " << name << std::endl;

	configuration.Set("data_directory", data_directory_);

	isDistributed_ = true;

	double error_variance_value_nodal;
	configuration.Set("error.variance_nodal", "v > 0", error_variance_value_nodal);
	std::cout << "error variance dist " << error_variance_value_nodal << std::endl;

	// Get pointer to the single instance of SimvascularGlobalArrayTransfer
	gat = SimvascularGlobalArrayTransfer::Get();

	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {
		int actualIdx = (gat->pointerMapInt_["local index of unique nodes"])[unitIdx];
		actualIdx--;
		for (int kk = 0; kk < 3; kk++) {
			int obsFuncVal = (gat->pointerMapInt_["observation function displacement"])[kk * conpar.nshg + actualIdx];
			if (obsFuncVal > 0) { // if we want to observe on this node
				state_obs_index_.push_back(aug_state_part.getDuplicatedStateIndex(ncounter));
				error_variance_value_.push_back(error_variance_value_nodal);
				icounter++;
			}
			ncounter++;
		}
	}

	Nobservation_local_ = icounter;

	load_data_from_datfile_ = 1;
	data_file_name_ = "saved_observations_local_nodaldisp";
	std::stringstream s_temp;
	s_temp << rank_;
	data_file_name_ = data_directory_ + "/" + data_file_name_ + s_temp.str() + ".dat";
	saved_obs_file_name_ = "saved_observations_local_nodaldisp" + s_temp.str() + ".dat";

	// some console output
	std::cout << "[" << name_ << "] at rank (" << rank_ << ")" << " # of obs.: " << Nobservation_local_ << std::endl;
}

void SimvascularNodalDisplacementObservation::ApplyOperator(const state& x, observation& Hx1, observation& Hx2, int obs_start_index) const {

	int icounter = 0;

	for (int kk = 0; kk < Nobservation_local_; kk++) {

		Hx1(obs_start_index+kk) = x(state_obs_index_[icounter]);
		Hx2(obs_start_index+kk) = Hx1(obs_start_index+kk);

		icounter++;
	}

}


/* */
/* */
SimvascularDistanceObservation::SimvascularDistanceObservation() {
}

SimvascularDistanceObservation::~SimvascularDistanceObservation() {
}

void SimvascularDistanceObservation::Initialize(std::string name, const SimvascularAugStatePart &aug_state_part, VerdandiOps &configuration) {

	int icounter = 0, ncounter = 0;

	name_ = name;
	std::cout << "Initializing " << name << std::endl;

	configuration.Set("data_directory", data_directory_);

	isDistributed_ = true;

	double error_variance_value_dist;
	configuration.Set("error.variance_dist", "v > 0", error_variance_value_dist);
	std::cout << "error variance p avg " << error_variance_value_dist << std::endl;

	// Get pointer to the single instance of SimvascularGlobalArrayTransfer
	gat = SimvascularGlobalArrayTransfer::Get();

	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

		int actualIdx = (gat->pointerMapInt_["local index of unique nodes"])[unitIdx];
		actualIdx--;

		int obsFuncVal = (gat->pointerMapInt_["observation function distance"])[actualIdx];

		if (obsFuncVal > 0) {
			//state_obs_index_.push_back(ncounter);
			error_variance_value_.push_back(error_variance_value_dist);
			icounter++;
		}

		ncounter++;
	}

	Nobservation_local_ = icounter;

	for (std::size_t jj=0; jj < aug_state_part.getSize(); jj++) {
		disp_duplicated_state_index_.push_back( aug_state_part.getDuplicatedStateIndex((int)jj) );
	}

	load_data_from_datfile_ = 0;


	//data_file_name_ = "saved_observations_local_walldistances";
	std::stringstream s_temp;
	s_temp << rank_;
	//data_file_name_ = data_directory_ + "/" + data_file_name_ + s_temp.str() + ".dat";

	saved_obs_file_name_ = "saved_observations_local_walldistances" + s_temp.str() + ".dat";

	// some console output
	std::cout << "[" << name_ << "] at rank (" << rank_ << ")" << " # of obs.: " << Nobservation_local_ << std::endl;
}

void SimvascularDistanceObservation::ApplyOperator(const state& x, observation& Hx1, observation& Hx2, int obs_start_index) const {

	int icounter = 0;

	// populate temporary array with displacement information from input state
	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

		int actualIdx = (gat->pointerMapInt_["local index of unique nodes"])[unitIdx];

		for(int varIdx=0; varIdx < 3; varIdx++) { // 3 components

			(gat->pointerMapDP_["temporary array"])[varIdx * conpar.nshg + actualIdx-1] = x(disp_duplicated_state_index_[icounter++]);

		}

	}

	// the values coming from the state vector need to communicated at the interprocessor boundaries
	estim_helpers_temp_comm();

	// call the distance computation function
	if (rank_ == 0)
		std::cout << "computing distance to wall data surfaces (obs. op.)" << std::endl;

	elmdist(gat->pointerMapDP_["temporary array"], // temporary array now has displacement from input state x
			gat->pointerMapDP_["coordinates"],
			gat->pointerMapDP_["distance"],
			gat->pointerMapDP_["distance normal"],
			gat->pointerMapDP_["M*distance"]);

	icounter = 0;

	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

		int actualIdx = (gat->pointerMapInt_["local index of unique nodes"])[unitIdx];
		actualIdx--;

		if ( gat->pointerMapInt_["observation function distance"][actualIdx] ) {
			Hx1(obs_start_index+icounter) = gat->pointerMapDP_["distance"][actualIdx];
			Hx2(obs_start_index+icounter) = gat->pointerMapDP_["M*distance"][actualIdx];
			icounter++;
		}

	}

}


/* */
/* */
SimvascularFlowPressObservation::SimvascularFlowPressObservation() {
}

SimvascularFlowPressObservation::~SimvascularFlowPressObservation() {
}

void SimvascularFlowPressObservation::Initialize(std::string name, const SimvascularAugStatePart &aug_state_part, VerdandiOps &configuration) {

	int icounter = 0, ncounter = 0;

	// Get pointer to the single instance of SimvascularGlobalArrayTransfer
	gat = SimvascularGlobalArrayTransfer::Get();

	name_ = name;
	std::cout << "Initializing " << name << std::endl;

	configuration.Set("data_directory", data_directory_);

	isDistributed_ = false;

	configuration.Set("Nobservation_flow","",0,Nobservation_flow_);
	configuration.Set("Nobservation_avgpressure","",0,Nobservation_avgpressure_);

	csobs_origins_.resize(Nobservation_flow_+Nobservation_avgpressure_);
	csobs_normals_.resize(Nobservation_flow_+Nobservation_avgpressure_);

	configuration.Set("csobs_origins",csobs_origins_);
	configuration.Set("csobs_normals",csobs_normals_);
	configuration.Set("csobs_radii",csobs_radii_);

	if ((unsigned int)Nobservation_flow_+Nobservation_avgpressure_ != csobs_origins_.size() ||
		(unsigned int)Nobservation_flow_+Nobservation_avgpressure_ != csobs_normals_.size() ||
		(unsigned int)Nobservation_flow_+Nobservation_avgpressure_ != csobs_radii_.size()) {

		throw Error("SimvascularObservationManager::Initialize: number of observation locations does not match listed locations!");
	}

	int error_variance_value_flow, error_variance_value_avgpress;

	configuration.Set("error.variance_flow", "v > 0", error_variance_value_flow);
	std::cout << "error.variance_flow " << error_variance_value_flow << std::endl;

	configuration.Set("error.variance_avgpress", "v > 0", error_variance_value_avgpress);
	std::cout << "error.variance_avgpress " << error_variance_value_avgpress << std::endl;

	if (rank_ == 0) {
		for (int kk = 0; kk < Nobservation_flow_+Nobservation_avgpressure_; kk++) {
			csobs_origins_[kk].Print();
			csobs_normals_[kk].Print();
			std::cout << csobs_radii_[kk] << std::endl;
		}
	}

	for (std::size_t jj=0; jj < aug_state_part.getSize(); jj++) {
		sol_duplicated_state_index_.push_back( aug_state_part.getDuplicatedStateIndex((int)jj) );
	}

	geom_points_ = vtkSmartPointer<vtkPoints>::New();
	geom_points_def_ = vtkSmartPointer<vtkPoints>::New();
	geom_ids_ = vtkSmartPointer<vtkIdList>::New();
	geom_UGrid_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
	geom_UGrid_def_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
	geom_surface_def_ = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();

	for (int unitIdx = 0; unitIdx < conpar.nshg; unitIdx++) {

		double coordVal[3];
		double dispVal[3];

		for (int kk = 0; kk < 3; kk++)
			coordVal[kk] = (gat->pointerMapDP_["coordinates"])[kk * conpar.nshg + unitIdx];

		for (int kk = 0; kk < 3; kk++)
			dispVal[kk] = (gat->pointerMapDP_["displacement"])[kk * conpar.nshg + unitIdx];

		geom_points_->InsertPoint(unitIdx,coordVal[0],coordVal[1],coordVal[2]);
		geom_points_def_->InsertPoint(unitIdx,coordVal[0]+dispVal[0],coordVal[1]+dispVal[1],coordVal[2]+dispVal[2]);

	}

	geom_UGrid_->SetPoints(geom_points_);
	geom_UGrid_def_->SetPoints(geom_points_def_);

	for (int kk = 0; kk < gat->global_mien.size(); kk++)
		for (int jj = 0; jj < gat->global_npro[kk]; jj++) {

			for (int ii = 0; ii < 4; ii++) {
				int nodeIndex = gat->global_mien[kk][ii * gat->global_npro[kk] + jj];
				geom_ids_->InsertNextId(--nodeIndex);
			}

			geom_UGrid_->InsertNextCell(VTK_TETRA,geom_ids_);

			geom_UGrid_def_->InsertNextCell(VTK_TETRA,geom_ids_);

			geom_ids_->Reset();
		}

	geom_surface_def_->SetInput(geom_UGrid_def_);

	// set up cross-sectional flow and pressure observation with VTK



	vtkSmartPointer<vtkDoubleArray> geom_vel_array = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> geom_pres_array = vtkSmartPointer<vtkDoubleArray>::New();

	geom_vel_array->SetNumberOfComponents(3);

	for (int unitIdx = 0; unitIdx < conpar.nshg; unitIdx++) {

		double solVal[4];

		for (int kk = 0; kk < 4; kk++)
			solVal[kk] = (gat->pointerMapDP_["solution"])[kk * conpar.nshg + unitIdx];

		geom_vel_array->InsertNextTuple3(solVal[0],solVal[1],solVal[2]);
		geom_pres_array->InsertNextTuple1(solVal[3]);
	}

	geom_UGrid_->GetPointData()->AddArray(geom_vel_array);
	geom_UGrid_->GetPointData()->GetArray(0)->SetName("velocity");
	geom_UGrid_->GetPointData()->AddArray(geom_pres_array);
	geom_UGrid_->GetPointData()->GetArray(1)->SetName("pressure");
	geom_UGrid_->Update();

	/* here we specify the number of cross-sections we want to observe flow on
	    	   the single slice plane may cut through several parts of the
	    	   mesh and only one part (i.e. the local geometry) is the desired cut.
	    	   the desired cut is not necssarily connected and may not be on the current processor

	    	   By default, vtkCutter generates triangulated surfaces
	 */

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

		/* experimental */
		geom_cutter_alt_ = vtkSmartPointer<vtkCutter>::New();
		geom_cutter_alt_->SetCutFunction(geom_plane_);
		geom_cutter_alt_->SetInput(geom_surface_def_->GetOutput());
		geom_cutter_alt_->Update();
		/*              */

		/*geom_connectivity_ = vtkSmartPointer<vtkConnectivityFilter>::New();
	    		geom_connectivity_->SetInputConnection(geom_cutter_->GetOutputPort());
	    		geom_connectivity_->SetExtractionModeToClosestPointRegion();
	    		geom_connectivity_->SetClosestPoint(csobs_origins_[kk](0),csobs_origins_[kk](1),csobs_origins_[kk](2));
	    		geom_connectivity_->Update();*/

		std::string fn1("obs_cut_"),fn2(".vtk");
		std::stringstream s_temp;

		s_temp << kk;
		fn1 = fn1 + s_temp.str() + fn2;

		geom_writer_ = vtkSmartPointer<vtkPolyDataWriter>::New();
		geom_writer_->SetFileName(fn1.c_str());
		geom_writer_->SetInput(geom_cutter_->GetOutput());
		geom_writer_->Write();

		/* experimental */
		fn1 = "obs_cut_surface_";
		fn1 = fn1 + s_temp.str() + fn2;
		geom_writer_->SetFileName(fn1.c_str());
		geom_writer_->SetInput(geom_cutter_alt_->GetOutput());
		geom_writer_->Write();
		/*              */

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
		//geom_connec_filters_.push_back(geom_connectivity_);

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
	}


	if (rank_ == numProcs_ - 1) { // note we want Nobservation_local_ to be nonzero on the last process
		for (int kk = 0; kk < Nobservation_flow_; ++kk) {
			Nobservation_local_++;
			error_variance_value_.push_back(error_variance_value_flow);
		}

		for (int kk = 0; kk < Nobservation_avgpressure_; ++kk) {
			Nobservation_local_++;
			error_variance_value_.push_back(error_variance_value_avgpress);
		}
	}

	load_data_from_datfile_ = 1;
	data_file_name_ = "saved_observations_single";
	data_file_name_ = data_directory_ + "/" + data_file_name_ + ".dat";
	saved_obs_file_name_ = "saved_observations_single.dat";

	// some console output
	std::cout << "[" << name_ << "] at rank (" << rank_ << ")" << " # of obs.: " << Nobservation_local_ << std::endl;

}

void SimvascularFlowPressObservation::ApplyOperator(const state& x, observation& Hx1, observation& Hx2, int obs_start_index) const {

	int icounter = 0;

	vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

	// Here we look through the state vector and grab the velocity values
	// to assign to a temporary array that will allow communication of those values to the
	// duplicated nodes on the interprocessor boundaries

	// the use of this temporary array could possibly
	// be avoided if model stateupdated is called
	// before the getinnovation function is called
	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

		int actualIdx = (gat->pointerMapInt_["local index of unique nodes"])[unitIdx];

		for(int varIdx=0; varIdx < 4; varIdx++) // ignore the 5th dof and beyond
			(gat->pointerMapDP_["temporary array"])[varIdx * conpar.nshg + actualIdx-1] = x(sol_duplicated_state_index_[icounter++]);;

	}

	// the values coming from the state vector need to communicated at the interprocessor boundaries
	estim_helpers_temp_comm();

	// now we reassign the velocities values to update the cross-sectional avg flow

	for (int unitIdx = 0; unitIdx < conpar.nshg; unitIdx++) {

		double solVal[4];

		for (int ii = 0; ii < 4; ii++)
			solVal[ii] = (gat->pointerMapDP_["temporary array"])[ii * conpar.nshg + unitIdx];

		geom_UGrid_->GetPointData()->GetArray(0)->SetTuple3(unitIdx,solVal[0],solVal[1],solVal[2]);
		geom_UGrid_->GetPointData()->GetArray(1)->SetTuple1(unitIdx,solVal[3]);

	}

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
		MPI_Reduce(&avgFlow_local, &avgFlow, 1, MPI_DOUBLE, MPI_SUM, numProcs_ - 1,iNewComm_C_);
		MPI_Reduce(&avgPres_local, &avgPres, 1, MPI_DOUBLE, MPI_SUM, numProcs_ - 1,iNewComm_C_);
		MPI_Reduce(   &area_local, &area,    1, MPI_DOUBLE, MPI_SUM, numProcs_ - 1,iNewComm_C_);

		if (rank_ == numProcs_ -1) { // note we only write to the observation on the last process
			(kk < Nobservation_flow_) ?
					Hx1(kk+obs_start_index) = avgFlow : Hx1(kk+obs_start_index) = avgPres / area;

			Hx2(kk+obs_start_index) = Hx1(kk+obs_start_index);
		}

	}

}


/////////////////////////////////
// CONSTRUCTORS AND DESTRUCTOR //
/////////////////////////////////


SimvascularObservationManager::SimvascularObservationManager()
:   Nobservation_(0),
    Nobservation_local_(0),
    Nobservation_nodal_(0),
    Nobservation_dist_(0),
    Nobservation_area_(0),
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


/*!
 \param[in] model model.
 \param[in] configuration_file configuration file.
 \tparam Model the model type; e.g. SimvascularVerdandiModel<double>
 */
template<class Model>
void SimvascularObservationManager::Initialize(const Model& model,
		string configuration_file) {
	VerdandiOps configuration(configuration_file);

	iNewComm_C_ = MPI_Comm_f2c(newcom.iNewComm);
	MPI_Comm_rank(iNewComm_C_, &rank_);
	MPI_Comm_size(iNewComm_C_, &numProcs_);

	Nstate_model_ = model.GetNstate();
	Nobservation_local_ = 0;

	// Get pointer to the single instance of SimvascularGlobalArrayTransfer
	gat = SimvascularGlobalArrayTransfer::Get();

	configuration.SetPrefix("observation.");
	//configuration.Set("use_restarts",use_restarts_);
	configuration.Set("data_directory", data_directory_);
	configuration.Set("Nskip", "v > 0", Nskip_);
	configuration.Set("initial_time", "", 0., initial_time_);
	configuration.Set("final_time", "", numeric_limits<double>::max(),final_time_);
	configuration.Set("data_period","",final_time_,data_period_);

	configuration.Set("execute_nodal_observations_","",0,execute_nodal_observations_);
	configuration.Set("execute_distance_observations_","",0,execute_distance_observations_);
	configuration.Set("Nobservation_flow","",0,Nobservation_flow_);
	configuration.Set("Nobservation_avgpressure","",0,Nobservation_avgpressure_);

	/* experimental */
	SimvascularObservationType* obsobjptr;

	if (execute_nodal_observations_ > 0 && nomodule.ideformwall > 0) {
		obsobjptr = new SimvascularNodalDisplacementObservation;
		obsobjptr->Initialize("Nodal",model.GetAugStateDstrb("displacement") , configuration);
		observations_dstrb_.push_back(obsobjptr);
	}

	if (execute_distance_observations_ && nomodule.ideformwall > 0) {
		obsobjptr = new SimvascularDistanceObservation;
		obsobjptr->Initialize("Distance",model.GetAugStateDstrb("displacement") , configuration);
		observations_dstrb_.push_back(obsobjptr);
	}

	if (Nobservation_flow_+Nobservation_avgpressure_ > 0) {
		obsobjptr = new SimvascularFlowPressObservation;
		obsobjptr->Initialize("FlowPress",model.GetAugStateDstrb("solution") , configuration);
		observations_single_.push_back(obsobjptr);
	}


	// count the number of observations on this processor
	for (auto it = observations_dstrb_.begin(); it != observations_dstrb_.end(); ++it)
		Nobservation_local_ += (*it)->getNobservation_local();

	for (auto it = observations_single_.begin(); it != observations_single_.end(); ++it)
		Nobservation_local_ += (*it)->getNobservation_local();

	// compute the global number of observation
    MPI_Allreduce(&Nobservation_local_, &Nobservation_, 1, MPI_INT, MPI_SUM, iNewComm_C_);

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

	int ncounter = 0;

	for (auto it = observations_dstrb_.begin(); it != observations_dstrb_.end(); ++it)
		for (int kk = 0; kk < (*it)->getNobservation_local(); ++kk) {
			error_variance_inverse_diag_.SetBuffer( local_start+ncounter, double(1) / (*it)->getErrorVarianceValue(kk)  );
			ncounter++;
		}

	if (rank_ == numProcs_ - 1) {

		for (auto it = observations_single_.begin(); it != observations_single_.end(); ++it)
			for (int kk = 0; kk < (*it)->getNobservation_local(); ++kk) {
				error_variance_inverse_diag_.SetBuffer( local_start+ncounter, double(1) / (*it)->getErrorVarianceValue(kk)  );
				ncounter++;
			}

	}

	error_variance_inverse_diag_.Flush();

#else

	throw ErrorUndefined("SimvascularObservationManager::"
	                                 "Initialize()", "Serial algorithm for storing observations not implemented.");

//#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
//	build_diagonal_sparse_matrix(Nobservation_, error_variance_value_,
//			error_variance_);
//	build_diagonal_sparse_matrix(Nobservation_,
//			double(double(1) / error_variance_value_),
//			error_variance_inverse_);
//#else
//	error_variance_.Reallocate(Nobservation_, Nobservation_);
//	error_variance_.SetIdentity();
//	Mlt(error_variance_value_, error_variance_);
//	error_variance_inverse_.Reallocate(Nobservation_, Nobservation_);
//	error_variance_inverse_.SetIdentity();
//	Mlt(double(double(1)/ error_variance_value_), error_variance_inverse_);
//#endif

#endif

}


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


/*!
      \param[in] discard_observation if set to true, each observation will be
      used at most one time.

      The discard observation flag currently has no use in SimvascularObservationManager.
 */
void SimvascularObservationManager::DiscardObservation(bool discard_observation)
{
	discard_observation_ = discard_observation;
}


void SimvascularObservationManager::GetObservation(
		SimvascularObservationManager::observation& observation) {
	throw ErrorUndefined(
			"void SimvascularObservationManager::GetObservation(const state& x,"
			"SimvascularObservationManager::observation& observation)");
}

////////////////
// INNOVATION //
////////////////


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
				std::cout << "loading data at time " << lower_bound << std::endl;

			LoadObservationSingleLocal(lower_bound, dataarrays_lower_);

			current_lower_bound_ = lower_bound;
		}

		if (upper_bound != current_upper_bound_) {

			if (rank_ == 0)
				std::cout << "loading data at time " << upper_bound << std::endl;

			LoadObservationSingleLocal(upper_bound, dataarrays_upper_);

			current_upper_bound_ = upper_bound;
		}
	}

	//dataarrays_lower_.Print();

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
    			t_alpha*dataarrays_lower_(kk) +
    			(1-t_alpha)*dataarrays_upper_(kk) ;

    	zHx2(kk) = -zHx2(kk) +
    			t_alpha*dataarrays_lower_(kk) +
    			(1-t_alpha)*dataarrays_upper_(kk) ;

    	//cout << "rank " << rank_ << " zhx1 " << zHx1(kk) << " zhx2 " << zHx2(kk) << endl;
    	//cout << "rank " << rank_ << " l " << dataarrays_lower_(kk) << " ";
    	//cout << "rank " << rank_ << " u " << dataarrays_upper_(kk) << endl;

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

    	//cout << "rank " << rank_ << "err.var.inv.diag.[" << kk << "]" << error_variance_inverse_diag_(local_start+kk) << endl;

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
    		iNewComm_C_);

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
                   iNewComm_C_);

    MPI_Allgatherv(zHx2.GetData(), Nobservation_local_, MPI_DOUBLE,
    		       innovation_fe.GetData(), obs_recvcount, obs_displ, MPI_DOUBLE,
    		       iNewComm_C_);

    delete [] obs_recvcount;
    delete [] obs_displ;

#endif
}


////////////
// ACCESS //
////////////


/*!
 \param[in] time a given time.
 */
bool SimvascularObservationManager::HasObservation(double time) {
	return time_ <= final_time_ && time_ >= initial_time_;
}


bool SimvascularObservationManager::HasObservation() const {
	return time_ <= final_time_ && time_ >= initial_time_;
}


/*!
 \return The total number of observation at current time.
 */
int SimvascularObservationManager::GetNobservation() const {
	return Nobservation_;
}


/*!
      \return The size of the local state vector.
 */
int SimvascularObservationManager::GetLocalNobservation() const {
	return Nobservation_local_; // this is used in reallocate routine for petsc matrices in ROUKF
}

///////////////
// OPERATORS //
///////////////


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


/*! This method is called after 'SetTime' set the time at which the
 operator is defined.
 \param[in] x a vector.
 \param[out] y the value of the operator applied to \a x. It is resized
 if needed.
 */
template<class state>
void SimvascularObservationManager::ApplyOperatorLocal(const state& x, observation& Hx1, observation& Hx2) {

	observation Hx_cs;

	//int icounter = 0;
	int ncounter = 0;

	int state_start, state_end;
	x.GetProcessorRange(state_start, state_end);

	Hx1.Reallocate(Nobservation_local_);
	Hx2.Reallocate(Nobservation_local_);

	Hx1.Zero();
	Hx2.Zero();

	for (auto it = observations_dstrb_.begin(); it != observations_dstrb_.end(); ++it) {
		(*it)->ApplyOperator(x,Hx1,Hx2,ncounter);
		ncounter += (*it)->getNobservation_local();
	}

	for (auto it = observations_single_.begin(); it != observations_single_.end(); ++it) {
		(*it)->ApplyOperator(x,Hx1,Hx2,ncounter);
		ncounter += (*it)->getNobservation_local();
	}

}


/*!
 \param[in] i row index.
 \param[in] j column index.
 \return The element (\a i, \a j) of the observation error variance.
 */
double SimvascularObservationManager::GetErrorVariance(int i, int j) const {
	throw ErrorUndefined("double SimvascularObservationManager"
			"::GetErrorVariance(int i, int j) const");
}


/*!
 \return The observation error covariance matrix.
 */
const SimvascularObservationManager::error_variance&
SimvascularObservationManager::GetErrorVariance() const {
	throw ErrorUndefined("const typename"
			"SimvascularObservationManager::error_variance& "
			"SimvascularObservationManager::GetErrorVariance() const");
}


/*!
 \return The inverse of the matrix of the observation error covariance.
 */
const SimvascularObservationManager::error_variance&
SimvascularObservationManager::GetErrorVarianceInverse() const {
	return error_variance_inverse_;
}


void SimvascularObservationManager::LoadObservationSingleLocal(int timeindex, Vector<double>& dataarray) {

	// in simvascular we assume that the time index is always an integral value
	// we need to sequentially read the file until
	int linetoread = timeindex / Nskip_;
	int ncounter = 0;

	for (auto it = observations_dstrb_.begin(); it != observations_dstrb_.end(); ++it) {
		(*it)->LoadData(linetoread,dataarray,ncounter);
		ncounter += (*it)->getNobservation_local();
	}

	if (rank_ == numProcs_ - 1)
		for (auto it = observations_single_.begin(); it != observations_single_.end(); ++it) {
			(*it)->LoadData(linetoread,dataarray,ncounter);
			ncounter += (*it)->getNobservation_local();
		}

}


/*! This function isn't actually used in any data assimilator routine

 */
template<class state>
void SimvascularObservationManager::SaveObservationSingleLocal(const state& x) {

	// we assumed that the time has been set from the model time
	if ((int)time_ % Nskip_ == 0) {

		//observation Hx1,Hx2;
		int ncounter = 0; //, ncounter2 = 0;

		if (rank_ == 0)
			std::cout << "SAVING OBSERVATIONS" << std::endl;

		for (auto it = observations_dstrb_.begin(); it != observations_dstrb_.end(); ++it) {
			(*it)->SaveObservations(x);
			ncounter += (*it)->getNobservation_local();
		}

		for (auto it = observations_single_.begin(); it != observations_single_.end(); ++it) {
			(*it)->SaveObservations(x);
			ncounter += (*it)->getNobservation_local();
		}

		if (rank_ == 0)
			std::cout << " [done] " << std::endl;

	}

}


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

#endif
