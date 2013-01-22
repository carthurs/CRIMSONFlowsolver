#include "GlobalArrayTransfer.h"

GlobalArrayTransfer* GlobalArrayTransfer::p_Instance_ = NULL;

GlobalArrayTransfer* GlobalArrayTransfer::Instance() {
	if (NULL == p_Instance_) {
		p_Instance_ = new GlobalArrayTransfer;
	}
	return p_Instance_;
}

GlobalArrayTransfer::GlobalArrayTransfer() {

}

GlobalArrayTransfer::~GlobalArrayTransfer() {

}

extern "C" void GlobalArrayAssignPointer(int* uniqptr, double* yoldptr, double* acoldptr, double* uoldptr,
		                                 double* coordptr, double* taptr,
		                                 double* distptr,
		                                 int* oyptr, int* oaptr, int* ouptr, int* odptr) {

	GlobalArrayTransfer::Instance()->global_inodesuniq_ptr = uniqptr;

	GlobalArrayTransfer::Instance()->global_yold_ptr = yoldptr;
	GlobalArrayTransfer::Instance()->global_acold_ptr = acoldptr;
	GlobalArrayTransfer::Instance()->global_uold_ptr = uoldptr;

	GlobalArrayTransfer::Instance()->global_coord_ptr = coordptr;
	GlobalArrayTransfer::Instance()->global_temporary_array_ptr = taptr;

	GlobalArrayTransfer::Instance()->global_xdist_ptr = distptr;

	GlobalArrayTransfer::Instance()->global_ilinobsfunc_sol_ptr = oyptr;
	GlobalArrayTransfer::Instance()->global_ilinobsfunc_acc_ptr = oaptr;
	GlobalArrayTransfer::Instance()->global_ilinobsfunc_disp_ptr = ouptr;
	GlobalArrayTransfer::Instance()->global_obsfunc_dist_ptr = odptr;

}

extern "C" void GlobalBlockedArrayAssignPointer(int npro_in, int nshl_in, int* ien_in) {
	GlobalArrayTransfer::Instance()->global_npro.push_back(npro_in);
	GlobalArrayTransfer::Instance()->global_nshl.push_back(nshl_in);
	GlobalArrayTransfer::Instance()->global_mien.push_back(ien_in);
}
