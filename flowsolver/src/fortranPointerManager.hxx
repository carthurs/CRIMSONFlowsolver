#ifndef _FORTRANPOINTERMANAGER_H_
#define _FORTRANPOINTERMANAGER_H_

#include <map>

class fortranBoundaryDataPointerManager {
public:

	static fortranBoundaryDataPointerManager* Get()
	{
		static fortranBoundaryDataPointerManager instance;
		return &instance;
	}

	std::map<int,double*> boundaryFlows;
	//std::map<int,double*> boundaryPressures;


};

#endif