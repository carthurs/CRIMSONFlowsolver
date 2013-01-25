#ifndef ROUKFMODIFIED_CXX

#include "ROUKFModified.hxx"
#include "Seldon.hxx"
#include "mpi.h"
using namespace Seldon;

template<class T, class Model, class ObservationManager>
ROUKFModified<T, Model, ObservationManager>::ROUKFModified() {

}

template<class T, class Model, class ObservationManager>
ROUKFModified<T, Model, ObservationManager>::~ROUKFModified() {

}

template<class T, class Model, class ObservationManager>
void ROUKFModified<T, Model, ObservationManager>::GetReducedStateErrorVariance(
		reduced_state_error_variance &Preduced) {

	Preduced.Reallocate(this->Nreduced_,this->Nreduced_);

	if (this->rank_ == this->model_.GetNumProcs() - 1) {

		reduced_state_error_variance L_temp, LU_inv_temp;

		L_temp.Reallocate(this->Nreduced_, this->Nreduced_);
		LU_inv_temp.Reallocate(this->Nreduced_, this->Nreduced_);

		L_temp.Fill(T(0));

		for (int i = 0; i < this->Nreduced_; i++)
			L_temp(i,i) = L_(i + this->Nstate_ - this->Nreduced_, i);

		MltAdd(T(1), SeldonNoTrans, L_temp, SeldonNoTrans, this->U_inv_, T(0), LU_inv_temp);

		MltAdd(T(1), SeldonNoTrans, LU_inv_temp, SeldonTrans, L_temp, T(0), Preduced);

		// broadcast contents of Preduced to all ranks

	}

	MPI_Bcast( Preduced.GetData(), Preduced.GetDataSize(), MPI_DOUBLE, this->model_.GetNumProcs()-1, MPI_COMM_WORLD );

}

#define ROUKFMODIFIED_CXX
#endif
