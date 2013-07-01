#ifndef _PHGLOBALARRAYTRANSFER_H_
#define _PHGLOBALARRAYTRANSFER_H_

#include <stdlib.h>
#include <vector>

using namespace std;

class PhGlobalArrayTransfer {
public:

	static PhGlobalArrayTransfer* Instance();
	~PhGlobalArrayTransfer();

	int* global_inodesuniq_ptr;

	double* global_yold_ptr;
	double* global_acold_ptr;
	double* global_uold_ptr;

	double* global_coord_ptr;
	double* global_temporary_array_ptr;

	double* global_xdist_ptr;
	double* global_df_fem_ptr;

	int* global_ilinobsfunc_sol_ptr;
	int* global_ilinobsfunc_acc_ptr;
	int* global_ilinobsfunc_disp_ptr;
	int* global_obsfunc_dist_ptr;

	double* global_lumped_parameter_P;
	double* global_lumped_parameter_Q;
	double* global_lumped_parameter_params;

	vector <int> global_npro;
	vector <int> global_nshl;
	vector <int*> global_mien;

	//int numblocks_;

private:
	PhGlobalArrayTransfer();

	static PhGlobalArrayTransfer* p_Instance_;
};

#endif
