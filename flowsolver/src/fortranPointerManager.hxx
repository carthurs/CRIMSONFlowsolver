#ifndef _FORTRANPOINTERMANAGER_H_
#define _FORTRANPOINTERMANAGER_H_

#include <map>
#include <cassert>
#include "gtest/gtest_prod.h"

class fortranBoundaryDataPointerManager {
	friend class testMultidom;
	friend class testMain;
	friend class testOrphans;
	friend class testMainWithZeroDDomain;
	FRIEND_TEST(testOrphans,checkNetlistDetectsBad3DInterfaceComponentOrientation);
	FRIEND_TEST(testOrphans,checkNetlistDetectsBadComponentAt3DInterface);
public:

	static fortranBoundaryDataPointerManager* Get()
	{
		static fortranBoundaryDataPointerManager instance;
		return &instance;
	}

	double* getBoundaryFlows(int surfaceIndex);
	double* getBoundaryPressures(int surfaceIndex);
	void setBoundaryFlows(int surfaceIndex, double* flowPointer);
	void setBoundaryPressures(int surfaceIndex, double* pressPointer);
	
private:
	// Make the constructor private; it's only ever called as a static method
	// via the public Get().
	fortranBoundaryDataPointerManager()
	{
		m_hasBoundaryPressures = false;
		m_hasBoundaryFlows = false;
	};

	// Ban (via making private) the copy constructor
	fortranBoundaryDataPointerManager(const fortranBoundaryDataPointerManager &old);

	// Ban (via making private) the assignment operator
	fortranBoundaryDataPointerManager &operator=(const fortranBoundaryDataPointerManager &old);

	// Ban (via making private) the destructor:
	~fortranBoundaryDataPointerManager(){};

	void tearDown()
    {
    	m_boundaryFlows.clear();
    	m_boundaryPressures.clear();
    }

    bool pointerNotDuplicated(double* pointerAboutToBeStored);

    bool m_hasBoundaryPressures;
    bool m_hasBoundaryFlows;

    std::map<int,double*> m_boundaryFlows;
	std::map<int,double*> m_boundaryPressures;
};

#endif