#ifndef _SIMVASCULARGLOBALARRAYTRANSFER_H_
#define _SIMVASCULARGLOBALARRAYTRANSFER_H_

#include <stdlib.h>
#include <vector>
#include <string>
#include <map>
#include <cassert>
#include <boost/shared_array.hpp>
#include <iostream>

class SimvascularGlobalArrayTransfer {
public:

	//! returns the only instance of this class
	static SimvascularGlobalArrayTransfer *Get()
	{
		if (!instance) {
			instance = new SimvascularGlobalArrayTransfer();
		}
		return instance;
		// static SimvascularGlobalArrayTransfer instance;
		// return &instance;
	}

	void insertDoublePointerPair(double* pointer, const std::string& fieldName);
	void initialiseForRCRFiltering(const int numberOfRCRSurfaces);
	void setPointerToWindkesselProximalResistance(double* pointerToProximalResistance, const int indexAmongstRCRs);
	void setPointerToWindkesselDistalResistance(double* pointerToDistalResistance, const int indexAmongstRCRs);
	void setPointerToWindkesselCompilance(double* pointerToCompliance, const int indexAmongstRCRs);
	void setPointerToRCRSurfacePressure(double* pointerToSurfacePressure, const int indexAmongstRCRs);
	void setSynchronisationDisabled(const std::string& keyName);

	void setPointerToFilteredNetlistParameter(double* pointerToFilteredParameter, const std::string parameterNameTag);
	const std::map<std::string, double*> getRawPointersToNetlistParameters() const;

	//! vector storing the block sizes for each element block
	std::vector <int> global_npro;

	//! vector storing the number of shape functions in each element block
	std::vector <int> global_nshl;

	//! vector storing pointers to the IEN array for each element block
	std::vector <int*> global_mien;

	//! map between string name and integer pointer
	std::map<std::string, int*> pointerMapInt_;


	// overloading wrapper so we can use char* if necessary:
	double getValueFromPointerMapDP(const char* keyName, const int dataLocation)
	{
		const std::string convertedChar(keyName);
		return getValueFromPointerMapDP(convertedChar, dataLocation);	
	}
	double getValueFromPointerMapDP(const std::string& keyName, const int dataLocation)
	{
		// If SimvascularGlobalArrayTransfer is managing a pointer to this value (i.e. it's from C++),
		// get the value from its true location in the class that owns it (e.g. compliance in the RCR class)	
		synchronisePointerMapDPArraysIfNecessary_get(keyName, dataLocation);
		// std::cout << "Return 1: " << keyName << " " << pointerMapDP_.at(keyName)[dataLocation] << std::endl;
		return pointerMapDP_.at(keyName)[dataLocation];
	}

	double* getRawPointerToSpecificValueRelatedToPointerMapDP(const char* keyName, const int dataLocation)
	{
		const std::string convertedChar(keyName);
		return getRawPointerToSpecificValueRelatedToPointerMapDP(convertedChar, dataLocation);
	}
	double* getRawPointerToSpecificValueRelatedToPointerMapDP(const std::string keyName, const int dataLocation)
	{
		double* ptrToReturn;
		bool needToReturnUnderlyingPointerToCppClassData = synchronisationEnabled_.at(keyName);
		if (needToReturnUnderlyingPointerToCppClassData)
		{
			ptrToReturn = actualDataPointerMap_.at(keyName).at(dataLocation);
			std::cout << "Return 2: " << keyName << " " << *ptrToReturn << std::endl;
		}
		else // else the pointer points to data which is held in a contiguous array (not a C++ class), so return that ptr:
		{
			ptrToReturn = &pointerMapDP_.at(keyName)[dataLocation];
		}
		return ptrToReturn;
	}

	// overloading wrapper so we can use char* if necessary:
	void setValueInPointerMapDP(const char* keyName, const int dataLocation, const double valueToSet)
	{
		const std::string convertedChar(keyName);
		setValueInPointerMapDP(convertedChar, dataLocation, valueToSet);
	}
	void setValueInPointerMapDP(const std::string& keyName, const int dataLocation, const double valueToSet)
	{
		pointerMapDP_.at(keyName)[dataLocation] = valueToSet;
		// If SimvascularGlobalArrayTransfer is managing a pointer to this value (i.e. it's from C++),
		// push the set value out to its true location in the class that owns it (e.g. compliance in the RCR class)	
		synchronisePointerMapDPArraysIfNecessary_set(keyName, dataLocation);
	}


	// overloading wrapper so we can use char* if necessary:
	double* getPointerToFirstArrayEntryFromPointerMapDP(const char* keyName)
	{
		const std::string convertedChar(keyName);
		return getPointerToFirstArrayEntryFromPointerMapDP(convertedChar);	
	}
	double* getPointerToFirstArrayEntryFromPointerMapDP(const std::string& keyName)
	{
		// synchronisePointerMapDPArraysIfNecessary_get(keyName);
		// WARNING: this array access is not synchronised, because some changes to it 
		// appear to happen in fortran, and so we can't manage synchronisation from here
		// currently. Instead, we ban /needing/ synchronisation here.
		assert(!synchronisationEnabled_.at(keyName));
		std::cout << "Return 3: " << keyName << " " << pointerMapDP_.at(keyName) << std::endl;
		return pointerMapDP_.at(keyName);
	}

	static void tearDown()
	{
		// global_npro.clear();
		// global_nshl.clear();
		// global_mien.clear();
		// pointerMapInt_.clear();
		// pointerMapDP_.clear();
		// actualDataPointerMap_.clear();
		// synchronisationEnabled_.clear();
		delete instance;
		instance = NULL;
	}

private:
	SimvascularGlobalArrayTransfer()
	{
		m_numberOfRCRSurfacesHasBeenSet = false;
	}
	~SimvascularGlobalArrayTransfer() {}
	SimvascularGlobalArrayTransfer(const SimvascularGlobalArrayTransfer &) { }
	SimvascularGlobalArrayTransfer &operator=(const SimvascularGlobalArrayTransfer &) { return *this; }

	void setupArraysForWindkesselFiltering();

	void synchronisePointerMapDPArraysIfNecessary_get(const std::string& keyName, const int dataLocation)
	{
		if (synchronisationEnabled_.at(keyName))
		{
			// synchronise to actual underlying data:
			pointerMapDP_.at(keyName)[dataLocation] = *(actualDataPointerMap_.at(keyName).at(dataLocation));
		}
	}

	void synchronisePointerMapDPArraysIfNecessary_set(const std::string& keyName, const int dataLocation)
	{
		if (synchronisationEnabled_.at(keyName))
		{
			// synchronise to actual underlying data:
			 *(actualDataPointerMap_.at(keyName).at(dataLocation)) = pointerMapDP_.at(keyName)[dataLocation];
		}
	}
	boost::shared_array<double> allRCRParameters_;
	boost::shared_array<double> interfacePressures_;

	int m_numberOfRCRSurfaces;
	bool m_numberOfRCRSurfacesHasBeenSet;

	//! map between string name and double pointer
	std::map<std::string, double*> pointerMapDP_;

	// This is for synchronising the pointerMapDP_ with the actual raw
	// data in the case where pointerMapDP_ has been built using the C++
	// boundary conditions, which do not assemble and lay out their
	// parameters in one contiguous array, as the Kalman filter interface
	// expects.
	//
	// To get around this, whenever we need to access pointerMapDP_, we first
	// synchronise all its data with the actual boundary conditions, using the
	// pointers stored in actualDataPointerMap_.
	//
	// This is untidy, but allows us to interface with the existing code.
	// \todo improve this by getting rid of the need for synchronisation,
	// by working always with pointers, not with an array of data as in pointerMapDP_.
	// This may not be possible without making a much worse mess unless we
	// get rid of the Fortran boundary conditions entriely.
	std::map<std::string, std::vector<double*>> actualDataPointerMap_;
	std::map<std::string, bool> synchronisationEnabled_;

	// For the netlists, we can be tidier and work with just one array of pointers
	// to the actual data:
	std::map<std::string, double*> netlistActualDataPointerMap_;

	static SimvascularGlobalArrayTransfer* instance;
};

#endif
