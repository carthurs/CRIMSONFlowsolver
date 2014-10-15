#include "fortranPointerManager.hxx"
#include <iostream>


extern "C" void giveflowpointertocpp(int& surfaceIndex,double*& flowPointer) {
	fortranBoundaryDataPointerManager::Get()->boundaryFlows.insert(std::pair<int,double*>(surfaceIndex,flowPointer));
	std::cout<< "WRITING... got: " << surfaceIndex << " " << *flowPointer << std::endl;
}