#ifndef SIMVASCULAROBSERVATIONMANAGER_CXX

#include "SimvascularObservationManager.hxx"

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

	delete [] dataarrays_lower_;
	delete [] dataarrays_upper_;
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

	configuration.SetPrefix("observation.");
	configuration.Set("use_synthetic_data",synthetic_data_);
	configuration.Set("data_directory", data_directory_);
	configuration.Set("Nskip", "v > 0", Nskip_);
	configuration.Set("initial_time", "", 0., initial_time_);
	configuration.Set("final_time", "", numeric_limits<double>::max(),
			final_time_);
	configuration.Set("error.variance", "v > 0", error_variance_value_);

	isize_solution_ = conpar.nshg * 4; // ignore last dof

	if (nomodule.ideformwall > 0)
		isize_displacement_ = conpar.nshg * NSD ;
	else
		isize_displacement_ = 0;

	isize_nshg_ = conpar.nshg;

    //
    // allocate space for the data (right now it is the size of the full state)
    // these are distributed arrays
    // isize_solution_ and isize_displacement_ should be the same as the array sizes in restart,
    // or else we have had a problem
    //
    dataarrays_lower_ = new double[2*isize_solution_ + isize_displacement_];
    dataarrays_upper_ = new double[2*isize_solution_ + isize_displacement_];

    for (int kk = 0; kk < 2*isize_solution_ + isize_displacement_; kk++) {
    	dataarrays_lower_[kk] = 0.0;
    	dataarrays_upper_[kk] = 0.0;
    }

    soln_lower_ = &dataarrays_lower_[0];
    soln_upper_ = &dataarrays_upper_[0];

    acc_lower_ = &dataarrays_lower_[isize_solution_];
    acc_upper_ = &dataarrays_upper_[isize_solution_];

    if (nomodule.ideformwall > 0) {
    	disp_lower_ = &dataarrays_lower_[2*isize_solution_];
    	disp_upper_ = &dataarrays_upper_[2*isize_solution_];
    }

    //
    // do some preprocessing on the 'simple' observation operator
    // tally number of observations on local processor
    //
    int obsCounter = 0;
    int actualIdx;
    int actualnshg = conpar.nshguniq;
    int obsFuncVal;

    for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

    	//phS->GetValue(*phS->GetRequiredField("local index of unique nodes"), unitIdx, 0, actualIdx);
    	actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];
    	actualIdx--;

    	for (int kk = 0; kk < 4; kk++) {

    		//phS->GetValue(*phS->GetRequiredField("observation function solution"), actualIdx, kk, obsFuncVal);
    		obsFuncVal = (gat->global_ilinobsfunc_sol_ptr)[kk * conpar.nshg + actualIdx];

    		//if (linobs_soln_[kk*isize_nshg_+actualIdx] > 0) {
    		if (obsFuncVal > 0) {
    			StateObsIndex_.PushBack(kk+4*unitIdx);
    			DataArraysObsIndex_.PushBack(kk*isize_nshg_+actualIdx);
    			obsCounter++;
    		}
    	}
    }

    for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

    	//phS->GetValue(*phS->GetRequiredField("local index of unique nodes"), unitIdx, 0, actualIdx);
    	actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];
    	actualIdx--;

    	for (int kk = 0; kk < 4; kk++) {

    		//phS->GetValue(*phS->GetRequiredField("observation function time derivative of solution"), actualIdx, kk, obsFuncVal);
    		obsFuncVal = (gat->global_ilinobsfunc_acc_ptr)[kk * conpar.nshg + actualIdx];

    		//if (linobs_acc_[kk*isize_nshg_+actualIdx] > 0) {
    		if (obsFuncVal > 0) {
    			StateObsIndex_.PushBack(kk+4*unitIdx + actualnshg*4);
    			DataArraysObsIndex_.PushBack(kk*isize_nshg_+actualIdx + isize_solution_);
    			obsCounter++;
    		}
    	}
    }

    if (nomodule.ideformwall > 0) {

    	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

    		//phS->GetValue(*phS->GetRequiredField("local index of unique nodes"), unitIdx, 0, actualIdx);
    		actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];
    		actualIdx--;

    		for (int kk = 0; kk < 3; kk++) {

    			//phS->GetValue(*phS->GetRequiredField("observation function displacement"), actualIdx, kk, obsFuncVal);
    			obsFuncVal = (gat->global_ilinobsfunc_disp_ptr)[kk * conpar.nshg + actualIdx];

    			//if (linobs_disp_[kk*isize_nshg_+actualIdx] > 0) {
    			if (obsFuncVal > 0) {
    				StateObsIndex_.PushBack(kk+3*unitIdx + actualnshg*4 + actualnshg*4);
    				DataArraysObsIndex_.PushBack(kk*isize_nshg_+actualIdx + isize_solution_*2);
    				obsCounter++;
    			}
    		}
    	}

    }

    //Nobservation_ = obsCounter;
    Nobservation_local_ = obsCounter;

    //
    // compute the global number of observations
    //
    MPI_Allreduce(&Nobservation_local_, &Nobservation_, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    //
    // initial observation error variance matrix
    //

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

	// set up the vtk data structures
	// this is here only for test purposes
	geom_points_ = vtkSmartPointer<vtkPoints>::New();
	geom_ids_ = vtkSmartPointer<vtkIdList>::New();
	geom_UGrid_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
	geom_vel_array_ = vtkSmartPointer<vtkDoubleArray>::New();

	geom_vel_array_->SetNumberOfComponents(3);

	double coordVal1,coordVal2,coordVal3;
	double solVal1,solVal2,solVal3;
	int nodeIndex;

	for (int unitIdx = 0; unitIdx < conpar.nshg; unitIdx++)
	{
		coordVal1 = (gat->global_coord_ptr)[0 * conpar.nshg + unitIdx];
		coordVal2 = (gat->global_coord_ptr)[1 * conpar.nshg + unitIdx];
		coordVal3 = (gat->global_coord_ptr)[2 * conpar.nshg + unitIdx];

		solVal1 = (gat->global_yold_ptr)[0 * conpar.nshg + unitIdx];
		solVal2 = (gat->global_yold_ptr)[1 * conpar.nshg + unitIdx];
		solVal3 = (gat->global_yold_ptr)[2 * conpar.nshg + unitIdx];

		geom_points_->InsertPoint(unitIdx,coordVal1,coordVal2,coordVal3);

		geom_vel_array_->InsertNextTuple3(solVal1,solVal2,solVal3);
	}

	geom_UGrid_->SetPoints(geom_points_);

	for (int kk = 0; kk < gat->global_mien.size(); kk++)
		for (int jj = 0; jj < gat->global_npro[kk]; jj++) {

//			phS->GetValueBlock(kk,jj,0,nodeIndex);
			nodeIndex = gat->global_mien[kk][0 * gat->global_npro[kk] + jj];
			geom_ids_->InsertNextId(--nodeIndex);

			nodeIndex = gat->global_mien[kk][1 * gat->global_npro[kk] + jj];
			geom_ids_->InsertNextId(--nodeIndex);

			nodeIndex = gat->global_mien[kk][2 * gat->global_npro[kk] + jj];
			geom_ids_->InsertNextId(--nodeIndex);

			nodeIndex = gat->global_mien[kk][3 * gat->global_npro[kk] + jj];
			geom_ids_->InsertNextId(--nodeIndex);

			geom_UGrid_->InsertNextCell(VTK_TETRA,geom_ids_);

			geom_ids_->Reset();
		}

	geom_UGrid_->GetPointData()->AddArray(geom_vel_array_);
	geom_UGrid_->GetPointData()->GetArray(0)->SetName("velocity");

	geom_UGrid_->Update();

	//	set up VTKcutters
	Nobservation_flow_ = 1;

	for (int kk = 0; kk < Nobservation_flow_; kk++) {

		// set the location of the cutting plane
		geom_plane_ = vtkSmartPointer<vtkPlane>::New();
		geom_plane_->SetOrigin(0,0,0);
		geom_plane_->SetNormal(0,0,1);
		geom_planes_.push_back(geom_plane_);

		geom_cutter_ = vtkSmartPointer<vtkCutter>::New();
		geom_cutter_->SetCutFunction(geom_plane_);
		geom_cutter_->SetInput(geom_UGrid_);
		geom_cutters_.push_back(geom_cutter_);

		geom_connectivity_ = vtkSmartPointer<vtkConnectivityFilter>::New();
		geom_connectivity_->SetInputConnection(geom_cutter_->GetOutputPort());
		geom_connectivity_->SetExtractionModeToClosestPointRegion();
		geom_connectivity_->SetClosestPoint(geom_plane_->GetOrigin());
		geom_connectivities_.push_back(geom_connectivity_);

	}

	flow_out_.open ("cross_section_mean_flow.dat");

//    vtkSmartPointer<vtkUnstructuredGridWriter> writer_ugrid = vtkUnstructuredGridWriter::New();
//    vtkSmartPointer<vtkPolyDataWriter> writer = vtkPolyDataWriter::New();
//    vtkSmartPointer<vtkUnstructuredGridWriter> writer_u = vtkUnstructuredGridWriter::New();
//    ostream *vtkout;
//
//    writer_ugrid->SetFileName("fromcons.vtk");
//    writer_ugrid->SetInput(tUGrid);
//    vtkout = writer_ugrid->OpenVTKFile();
//    writer_ugrid->Write();
//    writer_ugrid->CloseVTKFile(vtkout);
//
//    writer->SetFileName("fromcutter.vtk");
//    writer->SetInput(cutter->GetOutput());
//    vtkout = writer->OpenVTKFile();
//    writer->Write();
//    writer->CloseVTKFile(vtkout);
//
//    writer_u->SetFileName("fromconnectivity.vtk");
//    writer_u->SetInput(connectivityFilter->GetOutput());
//    vtkout = writer_u->OpenVTKFile();
//    writer_u->Write();
//    writer_u->CloseVTKFile(vtkout);

	// testing code for the CGAL stuff
//	int nodeIndex1,nodeIndex2,nodeIndex3,nodeIndex4;
//	double coordVal1,coordVal2,coordVal3;
//	const SCField *coord_field_ = phS->GetRequiredField("co-ordinates");
//
//	Point coord_pt;
//	vector<Point> coord_pts;
//	Polyhedron temp_poly;
//	Tree *temp_tree;
//
//	//cout << phS->GetNumBlocks() << endl;
//
//	for (int kk = 0; kk < coord_field_->GetNumUnits(); kk++) {
//		phS->GetValueDbl(*coord_field_,kk,0,coordVal1);
//		phS->GetValueDbl(*coord_field_,kk,1,coordVal2);
//		phS->GetValueDbl(*coord_field_,kk,2,coordVal3);
//		coord_pt = Point(coordVal1,coordVal2,coordVal3);
//		coord_pts.push_back(coord_pt);
//	}
//
//	for (int kk = 0; kk < phS->GetNumBlocks(); kk++) {
//		for (int jj = 0; jj < phS->GetBlockSize(kk); jj++) {
//
//			phS->GetValueIntBlock(kk,jj,0,nodeIndex1);
//			nodeIndex1--;
//
//            phS->GetValueIntBlock(kk,jj,1,nodeIndex2);
//            nodeIndex2--;
//
//            phS->GetValueIntBlock(kk,jj,2,nodeIndex3);
//            nodeIndex3--;
//
//            phS->GetValueIntBlock(kk,jj,3,nodeIndex4);
//            nodeIndex4--;
//
//            // add a tetrahedron
//            temp_poly.make_tetrahedron(coord_pts[nodeIndex1],coord_pts[nodeIndex2],coord_pts[nodeIndex3],coord_pts[nodeIndex4]);
//		}
//
//	}
//
//	// constructs AABB tree
//	temp_tree = new Tree(temp_poly.edges_begin(),temp_poly.edges_end());
//
//	// constructs plane query
//	CGAL_Vector vec(0.0,0.0,1.0);
//	Point a(0, 0, 0);
//	Plane plane_query(a,vec);
//
//	list<Object_and_primitive_id> intersections;
//
//	temp_tree->all_intersections(plane_query, back_inserter(intersections));
//
//	// gets intersection object
//	cout << "number of intersections " << intersections.size() << endl;
//
//	list<vector<double> > intersected_points;
//	vector<double> intersected_pt;
//	int pt_counter = 0;
//
//	list<Object_and_primitive_id>::iterator it;
//
//	for ( it=intersections.begin() ; it != intersections.end(); it++ ) {
//
//		CGAL::Object object = it->first;
//		Point testpoint;
//		if(CGAL::assign(testpoint,object)) {
////			cout << "intersection object is a point" << endl;
////			cout.precision(1);
////			cout << testpoint.x() << " " << testpoint.y() << " " << testpoint.z() << endl;
//
//			intersected_pt.push_back(testpoint.x());
//			intersected_pt.push_back(testpoint.y());
//			intersected_pt.push_back(testpoint.z());
//
//			intersected_points.push_back(intersected_pt);
//
//			intersected_pt.clear();
//
//            pt_counter++;
//		}
//	}
//
//	intersected_points.sort(compare_firstcoord);
//	intersected_points.unique(is_near());
//
//	intersected_points.sort(compare_secondcoord);
//	intersected_points.unique(is_near());
//
//	intersected_points.sort(compare_thirdcoord);
//	intersected_points.unique(is_near());
//
//	list<vector<double> >::iterator uniq_it;
//
//	for ( uniq_it=intersected_points.begin() ; uniq_it != intersected_points.end(); uniq_it++ ) {
//
//		intersected_pt = *uniq_it;
//		cout.precision(15);
//		cout << intersected_pt[0] << " " << intersected_pt[1] << " " << intersected_pt[2] << endl;
//	}
//
//	cout << "number of unique intersections " << intersected_points.size() << endl;
//
//	delete temp_tree;

	if (rank_ == 0)
		cout << "Simvascular Observation Manager initiated" << endl;

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

void SimvascularObservationManager::GetObservationFlow(
		SimvascularObservationManager::observation& observation) {


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
		if (rank_ == 0)
			cout << "loading data at time " << lower_bound << endl;

		if (synthetic_data_) {
			loadrestart(lower_bound,soln_lower_,acc_lower_,disp_lower_);
		}
		current_lower_bound_ = lower_bound;
	}

	if (upper_bound != current_upper_bound_) {
		if (upper_bound <= final_time_) {
			if (rank_ == 0)
				cout << "loading data at time " << upper_bound << endl;

			if (synthetic_data_) {
				loadrestart(upper_bound,soln_upper_,acc_upper_,disp_upper_);
			}
		}
		current_upper_bound_ = upper_bound;
	}

	//cout << "N observations " << Nobservation_ << endl;
	//cout << "N observations local " << Nobservation_local_ << endl;

	//
    // now we need to actually compute the innovation
    // compute the interpolation factors
	//

    double t_alpha = (time_ - current_lower_bound_)/(current_upper_bound_-current_lower_bound_);
    t_alpha = 1 - t_alpha;

    //
    // create petsc vector to scatter the innovation to all the processors.
    // for now, the innovation is required to be on
    //

    state zHx;
    int zHx_start, zHx_end;

    zHx.Reallocate(Nobservation_,Nobservation_local_);
    zHx.GetProcessorRange(zHx_start, zHx_end);

    //cout << "zHx_start " << zHx_start << endl;
    //cout << "zHx_end " << zHx_end << endl;

    // part of the innovation is on the local processor

    innovation.Reallocate(Nobservation_local_);

    int icounter = 0;

    int state_start, state_end;
    x.GetProcessorRange(state_start, state_end);

    for (int kk = zHx_start; kk < zHx_end; kk++) {

    	if (synthetic_data_) {
    		zHx.SetBuffer(kk, -x(state_start+StateObsIndex_(icounter)) +
    				t_alpha*dataarrays_lower_[DataArraysObsIndex_(icounter)] + (1-t_alpha)*dataarrays_upper_[DataArraysObsIndex_(icounter)] );
    	}

    	icounter++;
    }

    Vec zHx_seq;
    VecScatter ctx;

    VecScatterCreateToAll(zHx.GetPetscVector(), &ctx, &zHx_seq);
    VecScatterBegin(ctx, zHx.GetPetscVector(), zHx_seq,
    		INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx, zHx.GetPetscVector(), zHx_seq,
    		INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(&ctx);

    PetscInt ix[Nobservation_];
    PetscScalar zHx_seq_data[Nobservation_];
    Vector<int> ix_v;

    ix_v.SetData(Nobservation_, ix);
    ix_v.Fill();
    ix_v.Nullify();
    VecGetValues(zHx_seq, Nobservation_, ix, zHx_seq_data);
    Vector<double> zHx_seq_v;
    zHx_seq_v.SetData(Nobservation_, zHx_seq_data);
    innovation.Copy(zHx_seq_v);
    zHx_seq_v.Nullify();

    //innovation.Print();
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
void SimvascularObservationManager::ApplyOperator(const state& x, observation& y) const {
	throw ErrorUndefined("const typename"
			"void SimvascularObservationManager::ApplyOperator(const state& x, observation& y) const");

}

template<class state>
void SimvascularObservationManager::ApplyOperatorFlow(const state& x) {

	double solVal1,solVal2,solVal3;
	double vel1[3],vel2[3],vel3[3];

	double *tempcoord,*tempvel;
	double coord1[3],coord2[3],coord3[3];
	double A[3], B[3], C[3], triArea, avgFlow, tempL;

	int actualIdx;

	double val;

    int icounter = 0;

	vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

	for(int unitIdx=0; unitIdx < conpar.nshguniq; unitIdx++) {

		actualIdx = (gat->global_inodesuniq_ptr)[unitIdx];

		for(int varIdx=0; varIdx < 4; varIdx++) { // ignore the 5th dof and beyond

			val = x(icounter++);

			//phS->SetValue(*temp_field, actualIdx-1, varIdx, val);
			(gat->global_temporary_array_ptr)[varIdx * conpar.nshg + actualIdx-1] = val;

		}

	}

	// the values coming from the state vector need to communicated at the interprocessor boundaries
	if (numProcs_ > 1)
		temp_comm();

	for (int unitIdx = 0; unitIdx < conpar.nshg; unitIdx++)
	{
		solVal1 = (gat->global_temporary_array_ptr)[0 * conpar.nshg + unitIdx];
		solVal2 = (gat->global_temporary_array_ptr)[1 * conpar.nshg + unitIdx];
		solVal3 = (gat->global_temporary_array_ptr)[2 * conpar.nshg + unitIdx];

		geom_UGrid_->GetPointData()->GetArray("velocity")->SetTuple3(unitIdx,solVal1,solVal2,solVal3);
	}

	for (int kk = 0; kk < Nobservation_flow_ ; kk++) {
		geom_cutters_[kk]->Modified();
		geom_connectivities_[kk]->Update();


		//	for (int kk = 0; kk < geom_connectivity_->GetOutput()->GetNumberOfPoints(); kk++) {
		//		coord1 = geom_connectivity_->GetOutput()->GetPoint(kk);
		//		cout << kk << " " << coord1[0] << " " << coord1[1] << " " << coord1[2] << endl;
		//	}

		//	// loop through cells
		avgFlow = 0;
		for (int jj = 0; jj < geom_connectivities_[kk]->GetOutput()->GetNumberOfCells(); jj++) {
			geom_connectivities_[kk]->GetOutput()->GetCellPoints(jj,ptIds);

			tempvel = geom_connectivities_[kk]->GetOutput()->GetPointData()->GetArray("velocity")->GetTuple3(ptIds->GetId(0));
			vel1[0] = tempvel[0]; vel1[1] = tempvel[1]; vel1[2] = tempvel[2];

			tempvel = geom_connectivities_[kk]->GetOutput()->GetPointData()->GetArray("velocity")->GetTuple3(ptIds->GetId(1));
			vel2[0] = tempvel[0]; vel2[1] = tempvel[1]; vel2[2] = tempvel[2];

			tempvel = geom_connectivities_[kk]->GetOutput()->GetPointData()->GetArray("velocity")->GetTuple3(ptIds->GetId(2));
			vel3[0] = tempvel[0]; vel3[1] = tempvel[1]; vel3[2] = tempvel[2];

			tempcoord = geom_connectivities_[kk]->GetOutput()->GetPoint(ptIds->GetId(0));
			coord1[0] = tempcoord[0]; coord1[1] = tempcoord[1]; coord1[2] = tempcoord[2];

			tempcoord = geom_connectivities_[kk]->GetOutput()->GetPoint(ptIds->GetId(1));
			coord2[0] = tempcoord[0]; coord2[1] = tempcoord[1]; coord2[2] = tempcoord[2];

			tempcoord = geom_connectivities_[kk]->GetOutput()->GetPoint(ptIds->GetId(2));
			coord3[0] = tempcoord[0]; coord3[1] = tempcoord[1]; coord3[2] = tempcoord[2];

			//		cout << "cell " << kk << endl;
			//		cout << ptIds->GetId(0) << " " << ptIds->GetId(1) << " " << ptIds->GetId(2) << endl;
			//		cout << coord1[0] << " " << coord1[1] << " " << coord1[2] << endl;
			//		cout << coord2[0] << " " << coord2[1] << " " << coord2[2] << endl;
			//		cout << coord3[0] << " " << coord3[1] << " " << coord3[2] << endl;

			A[0] = coord2[0]-coord1[0];
			A[1] = coord2[1]-coord1[1];
			A[2] = coord2[2]-coord1[2];

			B[0] = coord3[0]-coord1[0];
			B[1] = coord3[1]-coord1[1];
			B[2] = coord3[2]-coord1[2];

			C[0] = (A[1]*B[2])-(B[1]*A[2]);
			C[1] = -(A[0]*B[2])+(B[0]*A[2]);
			C[2] = (A[0]*B[1])-(A[1]*B[0]);

			tempL = sqrt(C[0]*C[0]+C[1]*C[1]+C[2]*C[2]);

			triArea = 0.5*tempL;

			C[0] = C[0] / tempL;
			C[1] = C[1] / tempL;
			C[2] = C[2] / tempL;

			avgFlow += (double(1.0)/3.0)*(vel1[0]*C[0]+vel1[1]*C[1]+vel1[2]*C[2]+
					vel2[0]*C[0]+vel2[1]*C[1]+vel2[2]*C[2]+
					vel3[0]*C[0]+vel3[1]*C[1]+vel3[2]*C[2])*triArea;

		}

		if (rank_ == numProcs_ -1)
			this->flow_out_ << avgFlow << " ";

	}

	if (rank_ == numProcs_ -1)
		this->flow_out_ << endl;



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
	sprintf(restart_filename, "restart.%d.%d", timeindex, rank_+1);
	strcat(filename,"/");
	strcat(filename,restart_filename);

	//cout << "trying to load " << filename << endl;

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

	if (nomodule.ideformwall > 0) {
		readheader_(&irestart, "displacement?", (void*) iarray, &ithree, "double",
				iformat);
		isize = iarray[0]*iarray[1];
		readdatablock_(&irestart, "displacement?", (void*) disp, &isize,
				"double", iformat);
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
