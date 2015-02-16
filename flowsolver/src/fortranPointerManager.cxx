#include "fortranPointerManager.hxx"
#include <iostream>


extern "C" void giveflowpointertocpp(int& surfaceIndex,double*& flowPointer) {
	fortranBoundaryDataPointerManager::Get()->setBoundaryFlows(surfaceIndex,flowPointer);

}

extern "C" void givepressurepointertocpp(int& surfaceIndex,double*& pressPointer) {
	fortranBoundaryDataPointerManager::Get()->setBoundaryPressures(surfaceIndex,pressPointer);
}