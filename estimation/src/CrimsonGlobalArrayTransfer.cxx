#include "CrimsonGlobalArrayTransfer.h"
#include <iostream>
#include <cassert>
#include <boost/shared_array.hpp>

CrimsonGlobalArrayTransfer* CrimsonGlobalArrayTransfer::instance = NULL;

extern "C" void PhGlobalBlockedArrayAssignPointer(int npro_in, int nshl_in, int* ien_in) {
	CrimsonGlobalArrayTransfer::Get()->global_npro.push_back(npro_in);
	CrimsonGlobalArrayTransfer::Get()->global_nshl.push_back(nshl_in);
	CrimsonGlobalArrayTransfer::Get()->global_mien.push_back(ien_in);
}

extern "C" void PhAssignPointerInt(int* ptrInt, char* fieldName) {
	std::string tempName(fieldName);
	CrimsonGlobalArrayTransfer::Get()->pointerMapInt_.insert(std::pair<std::string, int*>(tempName, ptrInt));
}

extern "C" void PhAssignPointerDP(double* ptrDP, char* fieldName) {
	std::string tempName(fieldName);
	// std::cout << "data: " << CrimsonGlobalArrayTransfer::Get()->pointerMapDP_.at(tempName) << std::endl;
	
	CrimsonGlobalArrayTransfer::Get()->insertDoublePointerPair(ptrDP, tempName);
	// arrays passed from Fortran contain all the data in one contiguous memory block, so no
	// synchronisation is required to gather the data before accessing it:
	CrimsonGlobalArrayTransfer::Get()->setSynchronisationDisabled(tempName);
}

void CrimsonGlobalArrayTransfer::insertDoublePointerPair(double* pointer, const std::string& fieldName)
{
	pointerMapDP_.insert(std::make_pair(fieldName, pointer));
}

void CrimsonGlobalArrayTransfer::setSynchronisationDisabled(const std::string& keyName)
{
	synchronisationEnabled_.insert(std::make_pair(keyName, false));
}

void CrimsonGlobalArrayTransfer::initialiseForRCRFiltering(const int numberOfRCRSurfaces)
{
	m_numberOfRCRSurfaces = numberOfRCRSurfaces;
	m_numberOfRCRSurfacesHasBeenSet = true;
	setupArraysForWindkesselFiltering();
}

void CrimsonGlobalArrayTransfer::setupArraysForWindkesselFiltering()
{
	assert(m_numberOfRCRSurfacesHasBeenSet);
	{
		std::vector<double*> pointersToActualModelRCRParameters(3*m_numberOfRCRSurfaces, NULL);
		std::string valueName("WindkesselRCR_Params");
		actualDataPointerMap_.insert(std::make_pair(valueName, pointersToActualModelRCRParameters));
		synchronisationEnabled_.insert(std::make_pair(valueName, true));
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
		synchronisationEnabled_.insert(std::make_pair(valueName, true));
	}

	{
		interfacePressures_ = boost::shared_array<double> (new double[m_numberOfRCRSurfaces]());
		std::string valueName("WindkesselRCR_P");
		pointerMapDP_.insert(std::make_pair(valueName, interfacePressures_.get()));
	}
}

// To support arbitrary Netlist component filtering:
void CrimsonGlobalArrayTransfer::setPointerToFilteredNetlistParameter(double* pointerToFilteredParameter, const std::string parameterNameTag)
{
	netlistActualDataPointerMap_.insert(std::make_pair(parameterNameTag, pointerToFilteredParameter));
}

const std::map<std::string, double*> CrimsonGlobalArrayTransfer::getRawPointersToNetlistParameters() const
{
	return netlistActualDataPointerMap_;
}

void CrimsonGlobalArrayTransfer::setPointerToWindkesselProximalResistance(double* pointerToProximalResistance, const int indexAmongstRCRs)
{
	assert(m_numberOfRCRSurfacesHasBeenSet);
	pointerMapDP_.at("WindkesselRCR_Params")[indexAmongstRCRs*3 + 0] = *pointerToProximalResistance;
	actualDataPointerMap_.at("WindkesselRCR_Params")[indexAmongstRCRs*3 + 0] = pointerToProximalResistance;
}

void CrimsonGlobalArrayTransfer::setPointerToWindkesselDistalResistance(double* pointerToDistalResistance, const int indexAmongstRCRs)
{
	assert(m_numberOfRCRSurfacesHasBeenSet);
	pointerMapDP_.at("WindkesselRCR_Params")[indexAmongstRCRs*3 + 2] = *pointerToDistalResistance;
	actualDataPointerMap_.at("WindkesselRCR_Params")[indexAmongstRCRs*3 + 2] = pointerToDistalResistance;
}

void CrimsonGlobalArrayTransfer::setPointerToWindkesselCompilance(double* pointerToCompliance, const int indexAmongstRCRs)
{
	assert(m_numberOfRCRSurfacesHasBeenSet);
	pointerMapDP_.at("WindkesselRCR_Params")[indexAmongstRCRs*3 + 1] = *pointerToCompliance;
	actualDataPointerMap_.at("WindkesselRCR_Params")[indexAmongstRCRs*3 + 1] = pointerToCompliance;
}

void CrimsonGlobalArrayTransfer::setPointerToRCRSurfacePressure(double* pointerToSurfacePressure, const int indexAmongstRCRs)
{
	assert(m_numberOfRCRSurfacesHasBeenSet);
	pointerMapDP_.at("WindkesselRCR_P")[indexAmongstRCRs] = *pointerToSurfacePressure;
	actualDataPointerMap_.at("WindkesselRCR_P")[indexAmongstRCRs] = pointerToSurfacePressure;
}
