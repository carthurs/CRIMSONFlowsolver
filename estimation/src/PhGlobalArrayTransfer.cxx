#include "PhGlobalArrayTransfer.h"

/*
PhGlobalArrayTransfer* PhGlobalArrayTransfer::p_Instance_ = NULL;

PhGlobalArrayTransfer* PhGlobalArrayTransfer::Instance() {
	if (NULL == p_Instance_) {
		p_Instance_ = new PhGlobalArrayTransfer;
	}
	return p_Instance_;
}
*/

/*
PhGlobalArrayTransfer::PhGlobalArrayTransfer()
    :global_inodesuniq_ptr(NULL),
     global_yold_ptr(NULL),
     global_acold_ptr(NULL),
     global_uold_ptr(NULL),
     global_coord_ptr(NULL),
     global_temporary_array_ptr(NULL),
     global_xdist_ptr(NULL),
     global_df_fem_ptr(NULL),
     global_ilinobsfunc_sol_ptr(NULL),
     global_ilinobsfunc_acc_ptr(NULL),
     global_ilinobsfunc_disp_ptr(NULL),
     global_obsfunc_dist_ptr(NULL),
     global_lumped_parameter_P(NULL),
     global_lumped_parameter_Q(NULL),
     global_lumped_parameter_params(NULL),
     global_lumped_parameter_pout(NULL)
{

}

PhGlobalArrayTransfer::~PhGlobalArrayTransfer() {

}
*/

/* These functions are to be called from fortran subroutines */
/* Fortran interfaces are defined in common.f90 */

extern "C" void PhGlobalArrayAssignPointer(int* uniqptr, double* yoldptr, double* acoldptr, double* uoldptr,
		                                 double* coordptr, double* taptr,
		                                 double* distptr, double* df_femptr,
		                                 int* oyptr, int* oaptr, int* ouptr, int* odptr) {

	PhGlobalArrayTransfer::Get()->global_inodesuniq_ptr = uniqptr;

	PhGlobalArrayTransfer::Get()->global_yold_ptr = yoldptr;
	PhGlobalArrayTransfer::Get()->global_acold_ptr = acoldptr;
	PhGlobalArrayTransfer::Get()->global_uold_ptr = uoldptr;

	PhGlobalArrayTransfer::Get()->global_coord_ptr = coordptr;
	PhGlobalArrayTransfer::Get()->global_temporary_array_ptr = taptr;

	PhGlobalArrayTransfer::Get()->global_xdist_ptr = distptr;
	PhGlobalArrayTransfer::Get()->global_df_fem_ptr = df_femptr;

	PhGlobalArrayTransfer::Get()->global_ilinobsfunc_sol_ptr = oyptr;
	PhGlobalArrayTransfer::Get()->global_ilinobsfunc_acc_ptr = oaptr;
	PhGlobalArrayTransfer::Get()->global_ilinobsfunc_disp_ptr = ouptr;
	PhGlobalArrayTransfer::Get()->global_obsfunc_dist_ptr = odptr;

}

extern "C" void PhGlobalBlockedArrayAssignPointer(int npro_in, int nshl_in, int* ien_in) {
	PhGlobalArrayTransfer::Get()->global_npro.push_back(npro_in);
	PhGlobalArrayTransfer::Get()->global_nshl.push_back(nshl_in);
	PhGlobalArrayTransfer::Get()->global_mien.push_back(ien_in);
}

extern "C" void PhGlobalLumpedParameterArrayAssignPointer(double* p_ptr, double* q_ptr, double* param_ptr, double* pout_ptr) {
	PhGlobalArrayTransfer::Get()->global_lumped_parameter_P = p_ptr;
	PhGlobalArrayTransfer::Get()->global_lumped_parameter_Q = q_ptr;
	PhGlobalArrayTransfer::Get()->global_lumped_parameter_params = param_ptr;
	PhGlobalArrayTransfer::Get()->global_lumped_parameter_pout = pout_ptr;
}
