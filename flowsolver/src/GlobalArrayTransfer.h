#ifndef _GLOBALARRAYTRANSFER_H_
#define _GLOBALARRAYTRANSFER_H_

#include <stdlib.h>
#include <vector>

using namespace std;

class GlobalArrayTransfer {
public:

	static GlobalArrayTransfer* Instance();
	~GlobalArrayTransfer();

	int* global_inodesuniq_ptr;

	double* global_yold_ptr;
	double* global_acold_ptr;
	double* global_uold_ptr;

	double* global_coord_ptr;
	double* global_temporary_array_ptr;

	int* global_ilinobsfunc_sol_ptr;
	int* global_ilinobsfunc_acc_ptr;
	int* global_ilinobsfunc_disp_ptr;

	vector <int> global_npro;
	vector <int> global_nshl;
	vector <int*> global_mien;

	int numblocks_;

private:
	GlobalArrayTransfer();

	static GlobalArrayTransfer* p_Instance_;
};

#endif
