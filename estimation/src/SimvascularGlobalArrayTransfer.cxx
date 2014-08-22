#include "SimvascularGlobalArrayTransfer.h"
#include <iostream>

extern "C" void PhGlobalBlockedArrayAssignPointer(int npro_in, int nshl_in, int* ien_in) {
	SimvascularGlobalArrayTransfer::Get()->global_npro.push_back(npro_in);
	SimvascularGlobalArrayTransfer::Get()->global_nshl.push_back(nshl_in);
	SimvascularGlobalArrayTransfer::Get()->global_mien.push_back(ien_in);
}

extern "C" void PhAssignPointerInt(int* ptrInt, char* fieldName) {
	std::string tempName(fieldName);
	SimvascularGlobalArrayTransfer::Get()->pointerMapInt_.insert(std::pair<std::string, int*>(tempName, ptrInt));
}

extern "C" void PhAssignPointerDP(double* ptrDP, char* fieldName) {
	std::cout << *ptrDP << " ";

	std::string tempName(fieldName);
	SimvascularGlobalArrayTransfer::Get()->pointerMapDP_.insert(std::pair<std::string, double*>(tempName, ptrDP));

	std::cout << *SimvascularGlobalArrayTransfer::Get()->pointerMapDP_[tempName] << std::endl;
	
}
