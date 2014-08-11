#ifndef _SIMVASCULARGLOBALARRAYTRANSFER_H_
#define _SIMVASCULARGLOBALARRAYTRANSFER_H_

#include <stdlib.h>
#include <vector>
#include <string>
#include <map>

class SimvascularGlobalArrayTransfer {
public:

	static SimvascularGlobalArrayTransfer *Get()
	{
		static SimvascularGlobalArrayTransfer instance;
		return &instance;
	}

	~SimvascularGlobalArrayTransfer() {}

	std::vector <int> global_npro;
	std::vector <int> global_nshl;
	std::vector <int*> global_mien;

	std::map<std::string, int*> pointerMapInt_;
	std::map<std::string, double*> pointerMapDP_;

private:
	SimvascularGlobalArrayTransfer() {}
	SimvascularGlobalArrayTransfer(const SimvascularGlobalArrayTransfer &) { }
	SimvascularGlobalArrayTransfer &operator=(const SimvascularGlobalArrayTransfer &) { return *this; }

};

#endif
