#ifndef ROUKFMODIFIED_CXX
#define ROUKFMODIFIED_CXX

#include "ROUKFModified.hxx"
#include "Seldon.hxx"
#include "mpi.h"
using namespace Seldon;

template<class T, class Model, class ObservationManager>
ROUKFModified<T, Model, ObservationManager>::ROUKFModified() {

}

template<class T, class Model, class ObservationManager>
ROUKFModified<T, Model, ObservationManager>::~ROUKFModified() {
	if (this->model_.GetRank() == 0) {
		Poutfile_.close();
		Goutfile_.close();
	}
}

template<class T, class Model, class ObservationManager>
void ROUKFModified<T, Model, ObservationManager>::Initialize(string configuration_file,
		bool initialize_model,bool initialize_observation_manager) {

	Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>::
	Initialize(configuration_file,initialize_model,initialize_observation_manager);

	obsGram_.Reallocate(this->Nreduced_, this->Nreduced_);
	obsGram_.Zero();

	// set up file output
	std::string Pfilename = "Pred",filename_ext = ".dat";
	std::string Gfilename = "Gram";

	Pfilename = Pfilename+filename_ext;
	Gfilename = Gfilename+filename_ext;

	if (this->model_.GetRank() == this->model_.GetNumProcs() - 1) {
		Poutfile_.open(Pfilename.c_str(),ofstream::app);
		Goutfile_.open(Gfilename.c_str(),ofstream::app);

		Poutfile_ << std::scientific << std::setprecision( std::numeric_limits<double>::digits10 );
		Goutfile_ << std::scientific << std::setprecision( std::numeric_limits<double>::digits10 );
	}

	// load previous estimates if starting from previous estimation
	// this still needs to be finished


	// output for initial step if starting from scratch
	this->ComputeReducedStateErrorVariance();

	this->model_.WriteEstimates();

	if (this->model_.GetRank() == this->model_.GetNumProcs() - 1) {
		Poutfile_ << this->model_.GetTime() << endl;
		Preduced_.WriteText(Poutfile_);
		Goutfile_ << this->model_.GetTime() << endl;
		obsGram_.WriteText(Goutfile_);
	}

	this->observation_manager_.SetTime(this->model_,this->model_.GetTime());
	this->observation_manager_.SaveObservationSingleLocal(this->model_.GetState());

}

template<class T, class Model, class ObservationManager>
void ROUKFModified<T, Model, ObservationManager>::ComputeReducedStateErrorVariance() {

	if (this->rank_ == this->model_.GetNumProcs() - 1) {

		Preduced_.Reallocate(this->Nreduced_,this->Nreduced_);

		reduced_state_error_variance L_theta, LU_inv_temp;

		L_theta.Reallocate(this->Nreduced_, this->Nreduced_);
		LU_inv_temp.Reallocate(this->Nreduced_, this->Nreduced_);

		L_theta.Zero();

		int start_ind_local, end_ind_local;

		this->L_.GetProcessorRowRange(start_ind_local, end_ind_local);

		for (int i = 0; i < this->Nreduced_; i++)
			L_theta(i,i) = this->L_(start_ind_local + i + this->model_.GetLocalReducedStart(), i);

		MltAdd(T(1), SeldonNoTrans, L_theta, SeldonNoTrans, this->U_inv_, T(0), LU_inv_temp);

		MltAdd(T(1), SeldonNoTrans, LU_inv_temp, SeldonTrans, L_theta, T(0), Preduced_);

		//L_theta.Print();
		//this->U_inv_.Print();

	}

	// broadcast contents of Preduced to all ranks
	// reconsider whether this is necessary?
	//MPI_Bcast( Preduced.GetData(), Preduced.GetDataSize(), MPI_DOUBLE, this->model_.GetNumProcs()-1, MPI_COMM_WORLD );
}

//! Predict step
template <class T, class Model, class ObservationManager>
void ROUKFModified<T, Model, ObservationManager> ::Forward() {

	Verdandi::MessageHandler::Send(*this, "all", "::Forward begin");

	if (this->sigma_point_type_ == "simplex")
	{
#if defined(VERDANDI_WITH_MPI)
		throw Verdandi::ErrorUndefined("ROUKFModified::"
				"ForwardCustom()", "Parallel algorithm not"
				" implemented.");
#else
		model_state x;
		x.Copy(this->model_.GetState()); // <----------------------|

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
		{
			stringstream s_temp;
			s_temp << this->rank_;
			string outname = "x_before_sample-"+s_temp.str()+".dat";
			ofstream outfile(outname.c_str(), ofstream::app);
			outfile.precision(std::numeric_limits<double>::digits10);

			string outnameb = "x_before_sample-"+s_temp.str()+".bin";
			ofstream outfileb(outnameb.c_str(), ofstream::app);
			outfileb.precision(15);

			Vector <T> local_X;
			T* local_temp;
			int nxtemp = x.GetLocalM();
			VecGetArray(x.GetPetscVector(), &local_temp);  // get data pointer from petsc vector
			local_X.SetData(nxtemp,local_temp);  // set pointer for local vector

			local_X.WriteText(outfile);
			local_X.Write(outfileb);
			outfile << endl;
			local_X.Nullify();
		}
#endif

		/*** Sampling ***/

		//this->I_.Print();

		model_state_error_variance LC;
		sigma_point_matrix C;
		C.Copy(this->U_inv_);
		GetCholesky(C);

		LC.Copy(this->L_);
		MltAdd(T(1), this->L_, C, T(0), LC);

		// Computes X_n^{(i)+}.
		Reallocate(this->X_i_, this->Nstate_, this->Nsigma_point_, this->model_);
		for (int i = 0; i < this->Nsigma_point_; i++)
			SetCol(x, i, this->X_i_);

		MltAdd(T(1), LC, this->I_, T(1), this->X_i_);

		/*** Prediction ***/

		// Computes X_{n + 1}^-.
		x.Fill(T(0));
		model_state x_col;
		Reallocate(x_col, x.GetM(), this->model_);
		for (int i = 0; i < this->Nsigma_point_; i++)
		{
			GetCol(this->X_i_, i, x_col);

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
			{
				stringstream s_temp,s_temp2;
				s_temp << this->rank_;
				s_temp2 << i;
				string outname = "x_part_before_apply.p"+s_temp2.str()+"-"+s_temp.str()+".dat";
				ofstream outfile(outname.c_str(), ofstream::app);
				outfile.precision(std::numeric_limits<double>::digits10);

				string outnameb = "x_part_before_apply.p"+s_temp2.str()+"-"+s_temp.str()+".bin";
				ofstream outfileb(outnameb.c_str(), ofstream::app);
				outfileb.precision(15);

				Vector <T> local_X;
				T* local_temp;
				int nxtemp = x_col.GetLocalM();
				VecGetArray(x_col.GetPetscVector(), &local_temp);  // get data pointer from petsc vector
				local_X.SetData(nxtemp,local_temp);  // set pointer for local vector

				local_X.WriteText(outfile);
				local_X.Write(outfileb);
				outfile << endl;
				local_X.Nullify();
			}
#endif


			this->model_.ApplyOperator(x_col, i + 1 == this->Nsigma_point_, false); // <----------------------|

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
			{
				stringstream s_temp,s_temp2;
				s_temp << this->rank_;
				s_temp2 << i;
				string outname = "x_part_after_apply.p"+s_temp2.str()+"-"+s_temp.str()+".dat";
				ofstream outfile(outname.c_str(), ofstream::app);
				outfile.precision(std::numeric_limits<double>::digits10);

				string outnameb = "x_part_after_apply.p"+s_temp2.str()+"-"+s_temp.str()+".bin";
				ofstream outfileb(outnameb.c_str(), ofstream::app);
				outfileb.precision(15);

				Vector <T> local_X;
				T* local_temp;
				int nxtemp = x_col.GetLocalM();
				VecGetArray(x_col.GetPetscVector(), &local_temp);  // get data pointer from petsc vector
				local_X.SetData(nxtemp,local_temp);  // set pointer for local vector

				local_X.WriteText(outfile);
				local_X.Write(outfileb);
				outfile << endl;
				local_X.Nullify();
			}
#endif

			Add(T(this->alpha_), x_col, x);
			SetCol(x_col, i, this->X_i_);
		}

		//x_col.Print();
		//this->X_i_.Print();
		//int a,b,c,d;
		//this->X_i_.GetProcessorRowRange(a, b);
		//cout << "X_i range" << endl;
		//cout << a << " " << b << endl;

		/*** Resampling ***/

		if (this->with_resampling_)
		{
#if defined(VERDANDI_WITH_PETSC)
			throw Verdandi::ErrorUndefined("ROUKFModified::"
					"Forward()", "'resampling 'option "
					"not support yet.");
#endif
		}

        this->model_.GetState().Copy(x); // <----------------------|
        this->model_.StateUpdated(); // <----------------------|
#endif
	}
	else
	{
		throw Verdandi::ErrorUndefined("ROUKFModified::"
				"Forward()", "Parallel algorithm not "
				"implemented yet for the 'no"
				" simplex' cases.");

	}

	Verdandi::MessageHandler::Send(*this, "model", "forecast");
	Verdandi::MessageHandler::Send(*this, "observation_manager", "forecast");
	Verdandi::MessageHandler::Send(*this, "driver", "forecast");

	Verdandi::MessageHandler::Send(*this, "all", "::Forward end");


}

// Correction step
template<class T, class Model, class ObservationManager>
void ROUKFModified<T, Model, ObservationManager>::Analyze() {

    Verdandi::MessageHandler::Send(*this, "all", "::Analyze begin");

#if defined(VERDANDI_ROUKF_PARALLEL_INNOVATION)

#if defined(VERDANDI_WITH_MPI)
	throw Verdandi::ErrorUndefined("ROUKFModified::"
			"AnalyseCustomII()", "Parallel algorithm not"
			" implemented.");
#else

	this->observation_manager_.SetTime(this->model_, this->model_.GetTime() + 1);

	if (!this->observation_manager_.HasObservation())
	{
		Verdandi::MessageHandler::Send(*this, "all", "::Analyze end");
		return;
	}

	this->Nobservation_  = this->observation_manager_.GetNobservation();

	if (this->option_display_["show_time"] && this->rank_ == 0)
		cout << "Performing ROUKF to reach time step ["
		<< this->model_.GetTime() + 1 << "]..." << endl;

	if (this->sigma_point_type_ == "simplex")
	{
		// Computes L_{n + 1}.
		MltAdd(T(this->alpha_), this->X_i_, this->I_trans_, T(0), this->L_);

		// check with if the number of local observations is greater than 0
		int localHasObs = (this->observation_manager_.GetLocalNobservation() > 0);

		// Computes [HX_{n+1}^{*}].
		model_state_error_variance Z_i_p_orig;   // matrix of stored innovations (distributed)
		Z_i_p_orig.Reallocate( this->Nobservation_, this->Nsigma_point_,
				this->observation_manager_.GetLocalNobservation(), this->Nsigma_point_ );

		model_state_error_variance Z_i_p_fe;    // matrix of stored innovations (distributed)
		// with appropriate entries premultiplied by measurement error covariance inverse
		// and observed state error variance
		Z_i_p_fe.Reallocate( this->Nobservation_, this->Nsigma_point_,
				this->observation_manager_.GetLocalNobservation(), this->Nsigma_point_ );

		model_state x_col;   // temporary model state (distributed)
		Reallocate(x_col, this->Nstate_, this->model_);
		//this->model_.StateUpdated();    // <----------------------| update internal model state

		model_state z_col_p_orig;  // temporary innovation (distributed)
		model_state z_col_p_fe;

		model_state z_p_orig;   // stores weighted mean of innovation E(Zi) (distributed)
		z_p_orig.Reallocate(this->Nobservation_, this->observation_manager_.GetLocalNobservation() );
		z_p_orig.Zero();

		if (this->innovation_computation_)
			for (int i = 0; i < this->Nsigma_point_; i++)
			{
				GetCol(this->X_i_, i, x_col);  // retrieve stored model state (aka particle)
				this->observation_manager_.GetInnovation(x_col, z_col_p_orig, z_col_p_fe); // <----| obtain innovation

				Add(T(this->alpha_), z_col_p_orig, z_p_orig);   // weighted average

				SetCol(z_col_p_orig, i, Z_i_p_orig);   // store innovation
				SetCol(z_col_p_fe, i, Z_i_p_fe);
			}
		else
			throw Verdandi::ErrorUndefined("ROUKFModified::"
					"Analyze()", "Alternative to innovation not"
					" implemented.");


#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
		{
			stringstream s_temp;
			s_temp << this->rank_;
			string outname = "Z_i_p_orig-"+s_temp.str()+".dat";
			ofstream outfile(outname.c_str(), ofstream::app);
			outfile.precision(std::numeric_limits<double>::digits10);

			string outnameb = "Z_i_p_orig-"+s_temp.str()+".bin";
			ofstream outfileb(outnameb.c_str(), ofstream::app);
			outfileb.precision(15);

			Matrix<T, General, ColMajor, MallocAlloc<T> > local_A;
			int na = Z_i_p_orig.GetN();  // get number of columns
			T *local_a;
			MatGetArray(Z_i_p_orig.GetPetscMatrix(), &local_a);  // get data pointer from petsc
			int nlocal_A;
			int mlocal_A;
			MatGetLocalSize(Z_i_p_orig.GetPetscMatrix(), &mlocal_A, &nlocal_A);  // get number of local rows
			local_A.SetData(mlocal_A, na, local_a);  // set data to matrix

			local_A.WriteText(outfile);
			local_A.Write(outfileb);
			local_A.Nullify();
		}
#endif


		model_state_error_variance HL_p_orig;   // HL matrix "observed state sensitivity" (distributed)
		HL_p_orig.Reallocate(this->Nobservation_, this->Nreduced_,
				this->observation_manager_.GetLocalNobservation(), this->Nreduced_ );
		model_state_error_variance HL_p_fe;     // W^-1 * HL matrix "observed state sensitivity" (distributed)
		                                        // constructed from innovation vector premultiplied by mass matrix
		                                        // and observed state error variance
		HL_p_fe.Reallocate(this->Nobservation_, this->Nreduced_,
				this->observation_manager_.GetLocalNobservation(), this->Nreduced_ );

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
		{
			stringstream s_temp;
			s_temp << this->rank_;
			string outname = "I_trans-"+s_temp.str()+".dat";
			ofstream outfile(outname.c_str(), ofstream::app);
			outfile.precision(std::numeric_limits<double>::digits10);

			string outnameb = "I_trans-"+s_temp.str()+".bin";
			ofstream outfileb(outnameb.c_str(), ofstream::app);
			outfileb.precision(15);

			this->I_trans_.WriteText(outfile);
			this->I_trans_.Write(outfileb);
		}
#endif


		MltAdd(T(this->alpha_), Z_i_p_orig, this->I_trans_, T(0), HL_p_orig);
		MltAdd(T(this->alpha_), Z_i_p_fe, this->I_trans_, T(0), HL_p_fe);

		Matrix<T, General, ColMajor, MallocAlloc<T> > local_A; // matrix containing local entries for W^-1*HL_fe
		Matrix<T, General, ColMajor, MallocAlloc<T> > local_B; // matrix containing local entries for HL
		Matrix<T, General, ColMajor, MallocAlloc<T> > local_C(this->Nreduced_,this->Nreduced_); // partial HLtrans * W^-1 * HL
		Matrix<T, General, ColMajor, MallocAlloc<T> > global_C(this->Nreduced_,this->Nreduced_); // HLtrans * W^-1 * HL

		T *local_c, *global_c;

		local_C.Zero();
		global_C.Zero();

		local_c = local_C.GetData();  // get the pointer from the local part
		global_c = global_C.GetData();  // get the pointer from the global part

		// do the following only if number of "process-local observations" is greater than 0
		if (localHasObs) {

			// W^-1*HL_fe
			int na = HL_p_fe.GetN();  // get number of columns
			T *local_a;
			MatGetArray(HL_p_fe.GetPetscMatrix(), &local_a);  // get data pointer from petsc
			int nlocal_A;
			int mlocal_A;
			MatGetLocalSize(HL_p_fe.GetPetscMatrix(), &mlocal_A, &nlocal_A);  // get number of local rows
			local_A.SetData(mlocal_A, na, local_a);  // set data to matrix

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
			{
				stringstream s_temp;
				s_temp << this->rank_;
				string outname = "local_A-"+s_temp.str()+".dat";
				ofstream outfile(outname.c_str(), ofstream::app);
				outfile.precision(std::numeric_limits<double>::digits10);

				string outnameb = "local_A-"+s_temp.str()+".bin";
				ofstream outfileb(outnameb.c_str(), ofstream::app);
				outfileb.precision(15);

				local_A.WriteText(outfile);
				local_A.Write(outfileb);
			}
#endif

			// HL
			int nb = HL_p_orig.GetN();
			T *local_b;
			MatGetArray(HL_p_orig.GetPetscMatrix(), &local_b);
			int nlocal_B;
			int mlocal_B;
			MatGetLocalSize(HL_p_orig.GetPetscMatrix(), &mlocal_B, &nlocal_B);
			local_B.SetData(mlocal_B, nb, local_b);

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
			{
				stringstream s_temp;
				s_temp << this->rank_;
				string outname = "local_B-"+s_temp.str()+".dat";
				ofstream outfile(outname.c_str(), ofstream::app);
				outfile.precision(std::numeric_limits<double>::digits10);

				string outnameb = "local_B-"+s_temp.str()+".bin";
				ofstream outfileb(outnameb.c_str(), ofstream::app);
				outfileb.precision(15);

				local_B.WriteText(outfile);
				local_B.Write(outfileb);
			}
#endif

			// HLtrans * W^-1 * HL
			// do the multiply with the local
			MltAdd(T(1), SeldonTrans, local_A, SeldonNoTrans, local_B, T(0), local_C);

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
			{
				stringstream s_temp;
				s_temp << this->rank_;
				string outname = "local_C-"+s_temp.str()+".dat";
				ofstream outfile(outname.c_str(), ofstream::app);
				outfile.precision(std::numeric_limits<double>::digits10);

				string outnameb = "local_C-"+s_temp.str()+".bin";
				ofstream outfileb(outnameb.c_str(), ofstream::app);
				outfileb.precision(15);

				local_C.WriteText(outfile);
				local_C.Write(outfileb);
			}
#endif
		}

		// sum up HLtrans * W^-1 * HL across processes
		MPI_Allreduce(local_c, global_c, this->Nreduced_*this->Nreduced_, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

////		// allocate memory at root
//		double *rbuf;
//		if (this->rank_ == 0) {
//			rbuf = new double[this->Nreduced_*this->Nreduced_*this->model_.GetNumProcs()];
//		}
//
//		// gather to root
//		MPI_Gather( local_c, this->Nreduced_*this->Nreduced_, MPI_DOUBLE,
//				      rbuf, this->Nreduced_*this->Nreduced_, MPI_DOUBLE,
//				      0, MPI_COMM_WORLD );
//
//		// summation at root
//		if (this->rank_ == 0) {
//			for (int kk = 0; kk < this->model_.GetNumProcs(); kk++) {
//				for (int jj = 0; jj < this->Nreduced_*this->Nreduced_; jj++) {
//					global_c[jj] += rbuf[this->Nreduced_*this->Nreduced_*kk + jj];
//				}
//			}
//
//			delete [] rbuf;
//		}
//
//		// broadcast to all
//		MPI_Bcast( global_c, this->Nreduced_*this->Nreduced_, MPI_DOUBLE,
//				     0, MPI_COMM_WORLD );


		// compute the observability Gramian matrix
		if (this->rank_ == this->model_.GetNumProcs() - 1) {

			reduced_state_error_variance L_theta_inv, M1;

			L_theta_inv.Reallocate(this->Nreduced_, this->Nreduced_);
			L_theta_inv.Zero();
			M1.Reallocate(this->Nreduced_, this->Nreduced_);

			int start_ind_local, end_ind_local, state_reduced_start_local;
			this->L_.GetProcessorRowRange(start_ind_local, end_ind_local);
			state_reduced_start_local = this->model_.GetLocalReducedStart();

			for (int i = 0; i < this->Nreduced_; i++)
				L_theta_inv(i,i) = this->L_(start_ind_local + i + state_reduced_start_local, i);

			GetInverse(L_theta_inv);

			MltAdd(T(1), SeldonNoTrans, global_C, SeldonNoTrans, L_theta_inv, T(0), M1);
			MltAdd(T(1), SeldonTrans, L_theta_inv, SeldonNoTrans, M1, T(1), this->obsGram_);

			//MltAdd(T(1), SeldonNoTrans, LU_inv_temp, SeldonTrans, L_theta, T(0), Preduced_);

			//L_theta.Print();
			//this->U_inv_.Print();

		}

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
		{
			stringstream s_temp;
			s_temp << this->rank_;
			string outname = "global_C-"+s_temp.str()+".dat";
			ofstream outfile(outname.c_str(), ofstream::app);
			outfile.precision(std::numeric_limits<double>::digits10);

			string outnameb = "global_C-"+s_temp.str()+".bin";
			ofstream outfileb(outnameb.c_str(), ofstream::app);
			outfileb.precision(15);

			global_C.WriteText(outfile);
			global_C.Write(outfileb);
		}
#endif

		// U (U factor of state covariance matrix)
		this->U_inv_.SetIdentity(); //   initialize to identity

		Add(T(1), global_C, this->U_inv_); //   compute updated U

		GetInverse(this->U_inv_); //   compute updated U inverse

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
		{
			stringstream s_temp;
			s_temp << this->rank_;
			string outname = "U_inv-"+s_temp.str()+".dat";
			ofstream outfile(outname.c_str(), ofstream::app);
			outfile.precision(std::numeric_limits<double>::digits10);

			string outnameb = "U_inv-"+s_temp.str()+".bin";
			ofstream outfileb(outnameb.c_str(), ofstream::app);
			outfileb.precision(15);

			this->U_inv_.WriteText(outfile);
			this->U_inv_.Write(outfileb);
		}
#endif

		Vector <T> local_Z_avg;   //   local part of weighted mean E[Zi]

		// HLtrans * W^-1 * E[Zi]
		Vector <T> local_HLtransWinvZ(this->Nreduced_), global_HLtransWinvZ(this->Nreduced_);

		// Uinv * HLtrans * W^-1 * E[Zi]
		Vector <T> global_reduced_innovation(this->Nreduced_);

		global_HLtransWinvZ.Zero();  // initialize to zero
		local_HLtransWinvZ.Zero();

		T *local_d, *global_d;

		local_d = local_HLtransWinvZ.GetData(); // get the pointer from the local part
		global_d = global_HLtransWinvZ.GetData();  // get the pointer from the global vector

		// do the following only if number of "process-local observations" is greater than 0
		if (localHasObs) {

			// HLtrans * W^-1 * E[Zi]
			int nz = z_p_orig.GetLocalM();   // get number of rows
			T* local_z;
			VecGetArray(z_p_orig.GetPetscVector(), &local_z);  // get data pointer from petsc vector
			local_Z_avg.SetData(nz,local_z);  // set pointer for local vector

			MltAdd(T(1), SeldonTrans, local_A, local_Z_avg, T(0), local_HLtransWinvZ); // do the multiply locally

		}

		// sum up HLtrans * Winv * HL * E(Zi) across processes
		MPI_Allreduce(local_d, global_d, this->Nreduced_, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

//		// allocate memory at root
//		if (this->rank_ == 0) {
//			rbuf = new double[this->Nreduced_*this->model_.GetNumProcs()];
//		}
//
//		// gather to root
//		MPI_Gather( local_d, this->Nreduced_, MPI_DOUBLE,
//				      rbuf, this->Nreduced_, MPI_DOUBLE,
//				      0, MPI_COMM_WORLD );
//
//        // sum at root
//		if (this->rank_ == 0) {
//			for (int kk = 0; kk < this->model_.GetNumProcs(); kk++) {
//				for (int jj = 0; jj < this->Nreduced_; jj++) {
//					global_d[jj] += rbuf[this->Nreduced_*kk + jj];
//				}
//			}
//
//			delete [] rbuf;
//		}
//
//		// broadcast to all
//		MPI_Bcast( global_d, this->Nreduced_, MPI_DOUBLE,
//				     0, MPI_COMM_WORLD );

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
		{
			stringstream s_temp;
			s_temp << this->rank_;
			string outname = "local_Z-avg-"+s_temp.str()+".dat";
			ofstream outfile(outname.c_str(), ofstream::app);
			outfile.precision(std::numeric_limits<double>::digits10);

			string outnameb = "local_Z-avg-"+s_temp.str()+".bin";
			ofstream outfileb(outnameb.c_str(), ofstream::app);
			outfileb.precision(15);

			local_Z_avg.WriteText(outfile);
			local_Z_avg.Write(outfileb);
		}
#endif
#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
		{
			stringstream s_temp;
			s_temp << this->rank_;
			string outname = "global_D-"+s_temp.str()+".dat";
			ofstream outfile(outname.c_str(), ofstream::app);
			outfile.precision(std::numeric_limits<double>::digits10);

			string outnameb = "global_D-"+s_temp.str()+".bin";
			ofstream outfileb(outnameb.c_str(), ofstream::app);
			outfileb.precision(15);

			global_HLtransWinvZ.WriteText(outfile);
			global_HLtransWinvZ.Write(outfileb);
		}
#endif
#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
		{
			stringstream s_temp;
			s_temp << this->rank_;
			string outname = "local_D-"+s_temp.str()+".dat";
			ofstream outfile(outname.c_str(), ofstream::app);
			outfile.precision(std::numeric_limits<double>::digits10);

			string outnameb = "local_D-"+s_temp.str()+".bin";
			ofstream outfileb(outnameb.c_str(), ofstream::app);
			outfileb.precision(15);

			local_HLtransWinvZ.WriteText(outfile);
			local_HLtransWinvZ.Write(outfileb);
		}
#endif

		// Uinv * HLtrans * W^-1 * E[Zi]
		Mlt(T(-1), this->U_inv_, global_HLtransWinvZ, global_reduced_innovation);

		if (localHasObs) {
			// release pointers so that the temporary matrices can be destroyed
			// without deleting the data
			local_A.Nullify();   // release data to matrix to avoid destruction
			local_B.Nullify();
			local_Z_avg.Nullify();
		}

		// Updates.
		model_state& x =  this->model_.GetState();   // <----------------------|
		MltAdd(T(1), this->L_, global_reduced_innovation, T(1), x);
		this->model_.StateUpdated();   // <----------------------|

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
		{
			stringstream s_temp;
			s_temp << this->rank_;
			string outname = "x_updated-"+s_temp.str()+".dat";
			ofstream outfile(outname.c_str(), ofstream::app);
			outfile.precision(std::numeric_limits<double>::digits10);

			string outnameb = "x_updated-"+s_temp.str()+".bin";
			ofstream outfileb(outnameb.c_str(), ofstream::app);
			outfileb.precision(15);

			Vector <T> local_X;
			T* local_temp;
			int nxtemp = x.GetLocalM();
			VecGetArray(x.GetPetscVector(), &local_temp);  // get data pointer from petsc vector
			local_X.SetData(nxtemp,local_temp);  // set pointer for local vector

			local_X.WriteText(outfile);
			local_X.Write(outfileb);
            outfile << endl;

			local_X.Nullify();
		}
#endif

	}
	else
	{
		throw Verdandi::ErrorUndefined("ROUKFModified::"
				"AnalyseCustomII()", "Algorithm not"
				" implemented for the 'no"
				" simplex' cases.");
	}

	if (this->option_display_["show_time"] & this->rank_ == 0)
		cout << " done." << endl;

	Verdandi::MessageHandler::Send(*this, "model", "analysis");
	Verdandi::MessageHandler::Send(*this, "observation_manager", "analysis");
	Verdandi::MessageHandler::Send(*this, "driver", "analysis");

	Verdandi::MessageHandler::Send(*this, "all", "::AnalyzeCustomII end");

#endif

#else

#if defined(VERDANDI_WITH_MPI)
	throw Verdandi::ErrorUndefined("ROUKFModified::"
			"AnalyseCustom()", "Parallel algorithm not"
			" implemented.");
#else

	this->observation_manager_.SetTime(this->model_, this->model_.GetTime());

	if (!this->observation_manager_.HasObservation())
	{

		Verdandi::MessageHandler::Send(*this, "all", "::Analyze end");
		return;
	}

	this->Nobservation_  = this->observation_manager_.GetNobservation();

	if (this->option_display_["show_time"] && this->rank_ == 0)
		cout << "Performing ROUKF at time step ["
		<< this->model_.GetTime() << "]..." << endl;



	if (this->sigma_point_type_ == "simplex")
	{
		// Computes L_{n + 1}.
		MltAdd(T(this->alpha_), this->X_i_, this->I_trans_, T(0), this->L_);

		// Computes [HX_{n+1}^{*}].
		sigma_point_matrix Z_i_trans_orig(this->Nsigma_point_, this->Nobservation_);
		sigma_point_matrix Z_i_trans_fe(this->Nsigma_point_, this->Nobservation_);

		Z_i_trans_orig.Zero();
		Z_i_trans_fe.Zero();

		model_state x_col;
		Reallocate(x_col, this->Nstate_, this->model_);
		//this->model_.StateUpdated();    // <----------------------|

		observation z_col_orig(this->Nobservation_), z_col_fe(this->Nobservation_), z_orig(this->Nobservation_);
		z_orig.Zero();

		z_col_orig.Zero();
		z_col_fe.Zero();

		if (this->innovation_computation_)
			for (int i = 0; i < this->Nsigma_point_; i++)
			{
				GetCol(this->X_i_, i, x_col);
				this->observation_manager_.GetInnovation(x_col, z_col_orig, z_col_fe); // <----------------------|

				Add(T(this->alpha_), z_col_orig, z_orig); // weighted average
				SetRow(z_col_orig, i, Z_i_trans_orig);

				SetRow(z_col_fe, i, Z_i_trans_fe);

			}
		else
			throw Verdandi::ErrorUndefined("ROUKFModified::"
					"AnalyseCustom()", "Alternative to innovation not"
					" implemented.");

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
		{
			stringstream s_temp;
			s_temp << this->rank_;
			string outname = "Z_i_trans_orig-"+s_temp.str()+".dat";
			ofstream outfile(outname.c_str(), ofstream::app);
			outfile.precision(std::numeric_limits<double>::digits10);

			Z_i_trans_orig.WriteText(outfile);
		}
#endif

		sigma_point_matrix HL_trans_orig(this->Nreduced_, this->Nobservation_);
		sigma_point_matrix HL_trans_fe(this->Nreduced_, this->Nobservation_);

		MltAdd(T(this->alpha_), SeldonTrans, this->I_trans_, SeldonNoTrans, Z_i_trans_orig,
				T(0), HL_trans_orig);

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
		{
			stringstream s_temp;
			s_temp << this->rank_;
			string outname = "HL_trans_orig-"+s_temp.str()+".dat";
			ofstream outfile(outname.c_str(), ofstream::app);
			outfile.precision(std::numeric_limits<double>::digits10);

			HL_trans_orig.WriteText(outfile);
		}
#endif

		MltAdd(T(this->alpha_), SeldonTrans, this->I_trans_, SeldonNoTrans, Z_i_trans_fe,
				T(0), HL_trans_fe);

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
		{
			stringstream s_temp;
			s_temp << this->rank_;
			string outname = "HL_trans_fe-"+s_temp.str()+".dat";
			ofstream outfile(outname.c_str(), ofstream::app);
			outfile.precision(std::numeric_limits<double>::digits10);

			HL_trans_fe.WriteText(outfile);
		}
#endif

		sigma_point_matrix working_matrix_po(this->Nreduced_, this->Nobservation_), tmp;

		if (this->observation_error_variance_ == "matrix_inverse")
			Mlt(HL_trans_fe, this->observation_manager_.GetErrorVarianceInverse(),
					working_matrix_po);
		else
			throw Verdandi::ErrorUndefined("ROUKFModified::"
					"AnalyseCustom()", "Alternative to error variance inverse not"
					" implemented.");

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
		{
			stringstream s_temp;
			s_temp << this->rank_;
			string outname = "working_matrix_po-"+s_temp.str()+".dat";
			ofstream outfile(outname.c_str(), ofstream::app);
			outfile.precision(std::numeric_limits<double>::digits10);

			working_matrix_po.WriteText(outfile);
		}
#endif

		this->U_inv_.SetIdentity();

		MltAdd(T(1), SeldonNoTrans, working_matrix_po, SeldonTrans, HL_trans_orig, T(1), this->U_inv_);

		GetInverse(this->U_inv_);

#ifdef VERDANDI_ROUKF_DEBUG_OUTPUT
		{
			stringstream s_temp;
			s_temp << this->rank_;
			string outname = "U_inv-"+s_temp.str()+".dat";
			ofstream outfile(outname.c_str(), ofstream::app);
			outfile.precision(std::numeric_limits<double>::digits10);

			this->U_inv_.WriteText(outfile);
		}
#endif

		tmp.Reallocate(this->Nreduced_, this->Nobservation_);

		observation reduced_innovation(this->Nreduced_);

		MltAdd(T(1), this->U_inv_, working_matrix_po, T(0), tmp);

		if (this->innovation_computation_)
			MltAdd(T(-1), tmp, z_orig, T(0), reduced_innovation);
		else
			throw Verdandi::ErrorUndefined("ROUKFModified::"
					"AnalyseCustom()", "Alternative to innovation not"
					" implemented.");

		// Updates.
		model_state& x =  this->model_.GetState();
		MltAdd(T(1), this->L_, reduced_innovation, T(1), x);
		this->model_.StateUpdated(); // <----------------------|

	}
	else
	{
		throw Verdandi::ErrorUndefined("ROUKFModified::"
				"AnalyseCustom()", "Algorithm not"
				" implemented for the 'no"
				" simplex' cases.");
	}

	if (this->option_display_["show_time"] & this->rank_ == 0)
		cout << " done." << endl;

	Verdandi::MessageHandler::Send(*this, "model", "analysis");
	Verdandi::MessageHandler::Send(*this, "observation_manager", "analysis");
	Verdandi::MessageHandler::Send(*this, "driver", "analysis");

	Verdandi::MessageHandler::Send(*this, "all", "::AnalyzeCustom end");

#endif


#endif
}

template<class T, class Model, class ObservationManager>
void ROUKFModified<T, Model, ObservationManager>::FinalizeStep() {

	Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>::FinalizeStep();

	this->model_.WriteEstimates();

	this->ComputeReducedStateErrorVariance();

	if (this->model_.GetRank() == this->model_.GetNumProcs() - 1) {
		Poutfile_ << this->model_.GetTime() << endl;
		Preduced_.WriteText(Poutfile_);
		Goutfile_ << this->model_.GetTime() << endl;
		obsGram_.WriteText(Goutfile_);
	}

	this->observation_manager_.SetTime(this->model_,this->model_.GetTime());
	this->observation_manager_.SaveObservationSingleLocal(this->model_.GetState());

}

template<class T, class Model, class ObservationManager>
void ROUKFModified<T, Model, ObservationManager>::Finalize() {

	Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>::Finalize();

	// write out the estimation error factors

	Matrix<T, General, ColMajor, MallocAlloc<T> > local_L;

	// L matrix
	int n = this->L_.GetN();
	T *local_l;
	MatGetArray(this->L_.GetPetscMatrix(), &local_l);
	int nlocal_L;
	int mlocal_L;
	MatGetLocalSize(this->L_.GetPetscMatrix(), &mlocal_L, &nlocal_L);
	local_L.SetData(mlocal_L, n, local_l);

	{
		stringstream s_temp;
		s_temp << this->rank_;
		string outname = "L_local-final-"+s_temp.str()+".dat";
		ofstream outfile(outname.c_str(), ofstream::app);
		outfile.precision(std::numeric_limits<T>::digits10);
		local_L.WriteText(outfile);
	}

	local_L.Nullify();

	// U matrix
	if (this->model_.GetRank() == this->model_.GetNumProcs() - 1) {
		string outname = "U-final.dat";
		ofstream outfile(outname.c_str(), ofstream::app);
		outfile.precision(std::numeric_limits<T>::digits10);
		this->U_.WriteText(outfile);

	}


}

#endif
