#ifndef _FORTRANPOINTERMANAGER_H_
#define _FORTRANPOINTERMANAGER_H_

#include <map>
#include <cassert>
#include "gtest/gtest_prod.h"

class fortranBoundaryDataPointerManager {
	friend class testMultidom;
	friend class testMain;
	friend class testOrphans;
	FRIEND_TEST(testOrphans,checkNetlistDetectsBad3DInterfaceComponentOrientation);
	FRIEND_TEST(testOrphans,checkNetlistDetectsBadComponentAt3DInterface);
public:

	static fortranBoundaryDataPointerManager* Get()
	{
		static fortranBoundaryDataPointerManager instance;
		return &instance;
	}

	double* getBoundaryFlows(int surfaceIndex)
	{
		assert(hasBoundaryFlows);
		return boundaryFlows.at(surfaceIndex);
	}
	double* getBoundaryPressures(int surfaceIndex)
	{
		assert(hasBoundaryPressures);
		return boundaryPressures.at(surfaceIndex);
	}

	void setBoundaryFlows(int surfaceIndex, double* flowPointer)
	{
		boundaryFlows.insert(std::pair<int,double*>(surfaceIndex,flowPointer));
		hasBoundaryFlows = true;
	}
	void setBoundaryPressures(int surfaceIndex, double* pressPointer)
	{
		boundaryPressures.insert(std::pair<int,double*>(surfaceIndex,pressPointer));
		hasBoundaryPressures = true;
	}
private:
	// Make the constructor private; it's only ever called as a static method
	// via the public Get().
	fortranBoundaryDataPointerManager()
	{
		hasBoundaryPressures = false;
		hasBoundaryFlows = false;
	};

	// Ban (via making private) the copy constructor
	fortranBoundaryDataPointerManager(const fortranBoundaryDataPointerManager &old);

	// Ban (via making private) the assignment operator
	fortranBoundaryDataPointerManager &operator=(const fortranBoundaryDataPointerManager &old);

	// Ban (via making private) the destructor:
	~fortranBoundaryDataPointerManager(){};

	void tearDown()
    {
    	boundaryFlows.clear();
    	boundaryPressures.clear();
    }

    bool hasBoundaryPressures;
    bool hasBoundaryFlows;

    std::map<int,double*> boundaryFlows;
	std::map<int,double*> boundaryPressures;
};

#endif