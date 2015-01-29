#ifndef _FORTRANPOINTERMANAGER_H_
#define _FORTRANPOINTERMANAGER_H_

#include <map>

class fortranBoundaryDataPointerManager {
	friend class testMultidom;
	friend class testMain;
	friend class testOrphans;
public:

	static fortranBoundaryDataPointerManager* Get()
	{
		static fortranBoundaryDataPointerManager instance;
		return &instance;
	}

	std::map<int,double*> boundaryFlows;
	std::map<int,double*> boundaryPressures;
private:
	// Make the constructor private; it's only ever called as a static method
	// via the public Get().
	fortranBoundaryDataPointerManager(){};

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
};

#endif