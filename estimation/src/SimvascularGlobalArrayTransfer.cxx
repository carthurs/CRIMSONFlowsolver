#include "SimvascularGlobalArrayTransfer.h"
#include <iostream>
#include <cassert>
#include <boost/shared_array.hpp>

SimvascularGlobalArrayTransfer* SimvascularGlobalArrayTransfer::instance = NULL;

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
	std::string tempName(fieldName);
	// std::cout << "data: " << SimvascularGlobalArrayTransfer::Get()->pointerMapDP_.at(tempName) << std::endl;
	
	SimvascularGlobalArrayTransfer::Get()->insertDoublePointerPair(ptrDP, tempName);
	// arrays passed from Fortran contain all the data in one contiguous memory block, so no
	// synchronisation is required to gather the data before accessing it:
	SimvascularGlobalArrayTransfer::Get()->setSynchronisationDisabled(tempName);
}

void SimvascularGlobalArrayTransfer::insertDoublePointerPair(double* pointer, const std::string& fieldName)
{
	pointerMapDP_.insert(std::make_pair(fieldName, pointer));
}

void SimvascularGlobalArrayTransfer::setSynchronisationDisabled(const std::string& keyName)
{
	synchronisatoinEnabled_.insert(std::make_pair(keyName, false));
}

void SimvascularGlobalArrayTransfer::initialiseForRCRFiltering(const int numberOfRCRSurfaces)
{
	m_numberOfRCRSurfaces = numberOfRCRSurfaces;
	m_numberOfRCRSurfacesHasBeenSet = true;
	setupArraysForWindkesselFiltering();
}

void SimvascularGlobalArrayTransfer::setupArraysForWindkesselFiltering()
{
	assert(m_numberOfRCRSurfacesHasBeenSet);
	{
		std::vector<double*> pointersToActualModelRCRParameters(3*m_numberOfRCRSurfaces, NULL);
		std::string valueName("WindkesselRCR_Params");
		actualDataPointerMap_.insert(std::make_pair(valueName, pointersToActualModelRCRParameters));
		synchronisatoinEnabled_.insert(std::make_pair(valueName, true));
	}

	{
		allRCRParameters_ = boost::shared_array<double> (new double[3*m_numberOfRCRSurfaces]());
		std::string valueName("WindkesselRCR_Params");
		pointerMapDP_.insert(std::make_pair(valueName, allRCRParameters_.get()));
	}

	// boost::shared_array<double*> distalRCRPressures(new double[m_numberOfRCRSurfaces]());
	// pointerMapDP_.insert(std::pair<std::string, double*>("WindkesselRCR_Pdist", distalRCRPressures.get()));

	{	
		std::vector<double*> pointersToActualModelInterfacePressures(m_numberOfRCRSurfaces, NULL);
		std::string valueName("WindkesselRCR_P");
		actualDataPointerMap_.insert(std::make_pair(valueName, pointersToActualModelInterfacePressures));
		synchronisatoinEnabled_.insert(std::make_pair(valueName, true));
	}

	{
		interfacePressures_ = boost::shared_array<double> (new double[m_numberOfRCRSurfaces]());
		std::string valueName("WindkesselRCR_P");
		pointerMapDP_.insert(std::make_pair(valueName, interfacePressures_.get()));
	}
}

void SimvascularGlobalArrayTransfer::setPointerToWindkesselProximalResistance(double* pointerToProximalResistance, const int indexAmongstRCRs)
{
	assert(m_numberOfRCRSurfacesHasBeenSet);
	pointerMapDP_.at("WindkesselRCR_Params")[indexAmongstRCRs*3 + 0] = *pointerToProximalResistance;
	actualDataPointerMap_.at("WindkesselRCR_Params")[indexAmongstRCRs*3 + 0] = pointerToProximalResistance;
}

void SimvascularGlobalArrayTransfer::setPointerToWindkesselDistalResistance(double* pointerToDistalResistance, const int indexAmongstRCRs)
{
	assert(m_numberOfRCRSurfacesHasBeenSet);
	pointerMapDP_.at("WindkesselRCR_Params")[indexAmongstRCRs*3 + 2] = *pointerToDistalResistance;
	actualDataPointerMap_.at("WindkesselRCR_Params")[indexAmongstRCRs*3 + 2] = pointerToDistalResistance;
}

void SimvascularGlobalArrayTransfer::setPointerToWindkesselCompilance(double* pointerToCompliance, const int indexAmongstRCRs)
{
	assert(m_numberOfRCRSurfacesHasBeenSet);
	pointerMapDP_.at("WindkesselRCR_Params")[indexAmongstRCRs*3 + 1] = *pointerToCompliance;
	actualDataPointerMap_.at("WindkesselRCR_Params")[indexAmongstRCRs*3 + 1] = pointerToCompliance;
}

void SimvascularGlobalArrayTransfer::setPointerToRCRSurfacePressure(double* pointerToSurfacePressure, const int indexAmongstRCRs)
{
	assert(m_numberOfRCRSurfacesHasBeenSet);
	pointerMapDP_.at("WindkesselRCR_P")[indexAmongstRCRs] = *pointerToSurfacePressure;
	actualDataPointerMap_.at("WindkesselRCR_P")[indexAmongstRCRs] = pointerToSurfacePressure;
}
