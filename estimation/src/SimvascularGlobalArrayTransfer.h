#ifndef _SIMVASCULARGLOBALARRAYTRANSFER_H_
#define _SIMVASCULARGLOBALARRAYTRANSFER_H_

#include <stdlib.h>
#include <vector>
#include <string>
#include <map>

class SimvascularGlobalArrayTransfer {
public:

	//! returns the only instance of this class
	static SimvascularGlobalArrayTransfer *Get()
	{
		static SimvascularGlobalArrayTransfer instance;
		return &instance;
	}

	~SimvascularGlobalArrayTransfer() {}

	//! vector storing the block sizes for each element block
	std::vector <int> global_npro;

	//! vector storing the number of shape functions in each element block
	std::vector <int> global_nshl;

	//! vector storing pointers to the IEN array for each element block
	std::vector <int*> global_mien;

	//! map between string name and integer pointer
	std::map<std::string, int*> pointerMapInt_;

	//! map between string name and double pointer
	std::map<std::string, double*> pointerMapDP_;

private:
	SimvascularGlobalArrayTransfer() {}
	SimvascularGlobalArrayTransfer(const SimvascularGlobalArrayTransfer &) { }
	SimvascularGlobalArrayTransfer &operator=(const SimvascularGlobalArrayTransfer &) { return *this; }

};

#endif
