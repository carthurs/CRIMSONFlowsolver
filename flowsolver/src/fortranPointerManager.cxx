#include "fortranPointerManager.hxx"
#include <iostream>


extern "C" void giveflowpointertocpp(int& surfaceIndex,double*& flowPointer) {
	fortranBoundaryDataPointerManager::Get()->boundaryFlows.insert(std::pair<int,double*>(surfaceIndex,flowPointer));

}

extern "C" void givepressurepointertocpp(int& surfaceIndex,double*& pressPointer) {
	fortranBoundaryDataPointerManager::Get()->boundaryPressures.insert(std::pair<int,double*>(surfaceIndex,pressPointer));
}