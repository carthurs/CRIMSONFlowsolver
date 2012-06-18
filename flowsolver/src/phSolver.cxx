/****************************************************************************** 

 (c) 2009-2010 Scientific Computation Research Center,
 Rensselaer Polytechnic Institute. All rights reserved.

 The LICENSE-SCOREC file included with this distribution describes the terms
 of the SCOREC Non-Commercial License this program is distributed under.

 *******************************************************************************/

#include "phSolver.h"

namespace phSolverMessagePrinter {
void printDebugMsg(char* debugMsg);
void printErrorMsg(char* errorMsg);
void printWarningMsg(char* warningMsg);
}

//global static pointer used to ensure single instance of the class
phSolver* phSolver::p_Instance_ = NULL;

phSolver* phSolver::Instance() {
	if (NULL == p_Instance_) {
		p_Instance_ = new phSolver;
	}
	return p_Instance_;
}

phSolver::phSolver() {
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs_);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

	workfc.numpe = numProcs_;
	workfc.myrank = rank_;

	this->setCommunicator(MPI_COMM_WORLD);
}

void phSolver::setCommunicator(MPI_Comm comm) {

	this->iNewComm_C = comm;
	newcom.iNewComm = MPI_Comm_c2f(this->iNewComm_C); // modifies newcom in fortran common block
	MPI_Comm_size(this->iNewComm_C, &numProcs_);
	MPI_Comm_rank(this->iNewComm_C, &rank_);

	workfc.numpe = numProcs_;
	workfc.myrank = rank_;
}

MPI_Comm phSolver::getCommunicator() {
	return this->iNewComm_C;
}

void phSolver::readConfiguration() {
	this->input_fform(); // read solver.inp
}

int phSolver::readMeshAndSolution_fromFiles(char* pathToProcsCaseDir) {

	strncpy(pc_pathToProcsCaseDir_, pathToProcsCaseDir, 1024);

	char cwd[1024];
	char errorMessage[256];
	if (getcwd(cwd, sizeof(cwd)) == NULL) {
		sprintf(errorMessage, "getcwd() failed");
		phSolverMessagePrinter::printErrorMsg(errorMessage);
		return 1;
	}
	if (-1 == chdir(pc_pathToProcsCaseDir_)) {
		sprintf(errorMessage, "chdir(%s) failed", pc_pathToProcsCaseDir_);
		phSolverMessagePrinter::printErrorMsg(errorMessage);
		return 1;
	}
	/* in readnblk the required and provided fields lists are registered */
	//proces_(); // read geombc and restart from file to populate mesh, solution, and other field arrays
	input(&numProcs_, &rank_);
	proces();
	if (-1 == chdir(cwd)) {
		sprintf(errorMessage, "chdir(%s) failed", cwd);
		phSolverMessagePrinter::printErrorMsg(errorMessage);
		return 1;
	}

	return 0;
}
//
///*
// find the field and update its data - this function is meant to be called from common/readnblk.f and other locations where data is read from geombc and restart files
// */
//// see if it is possible to avoid using a global function here
//extern "C" void phSolverUpdateField(char* pc_name, char* pc_domainType,
//		char* pc_valType, int i_numUnits, int i_numVar, double *pd_dataArray,
//		int *pi_dataArray) {
//	SCField inField(pc_name, pc_domainType, i_numUnits, i_numVar, pc_valType);
//
////	char debugMessage[256];
////	sprintf(debugMessage, "pointer to field %s: dbl: %d int: %i",
////			inField.GetName(), pd_dataArray, pi_dataArray);
////	phSolverMessagePrinter::printDebugMsg(debugMessage);
//
//	void* pv_dataArray = NULL;
//	if (pd_dataArray != NULL && pi_dataArray == NULL) {
//		pv_dataArray = (void*) pd_dataArray;
//	} else if (pd_dataArray == NULL && pi_dataArray != NULL) {
//		pv_dataArray = (void*) pi_dataArray;
//	} else {
////		char warningMsg[256];
////		sprintf(warningMsg, "pointer to field %s was not updated",
////				inField.GetName());
////		phSolverMessagePrinter::printWarningMsg(warningMsg);
//	}
//
//	if (1 == phSolver::Instance()->UpdateField(inField, pv_dataArray)) {
////		char warningMsg[256];
////		sprintf(warningMsg, "field %s was not updated", inField.GetName());
////		phSolverMessagePrinter::printWarningMsg(warningMsg);
//	}
//
//}

extern "C" void phSolverUpdateField(char* pc_name, char* pc_domainType,
		char* pc_valType, int i_numUnits, int i_numVar, double *pd_dataArray,
		int *pi_dataArray) {

}

int phSolver::SolverInit() {
	int myrank;
	MPI_Comm_rank(this->iNewComm_C, &myrank);

	char errorMsg[256], warningMsg[256];

	if (-1 == chdir(pc_pathToProcsCaseDir_)) {
		sprintf(errorMsg, "chdir(%s) failed", pc_pathToProcsCaseDir_);
		phSolverMessagePrinter::printErrorMsg(errorMsg);
		return 1;
	}
//	phSolverMessagePrinter::printDebugMsg("Initializing solver...");
	itrdrv_init(); // calls the fortran routine

	return 0;
}

void phSolver::SolverForwardInit() {
	itrdrv_iter_init(); // calls the fortran routine
}
void phSolver::SolverForwardStep() {
	itrdrv_iter_step(); // calls the fortran routine
}
void phSolver::SolverForwardFinalize() {
	itrdrv_iter_finalize(); // calls the fortran routine

}
void phSolver::SolverFinalize() {
	itrdrv_finalize(); // calls the fortran routine
}

//int phSolver::SetValue(SCField field, int unitIndex, int varIndex,
//		double valueIn) {
//	//printField(field);
//	double *dataArray = (double*) UMap_FieldToPtr_[field];
//	if (NULL != dataArray) {
//		dataArray[varIndex * field.GetNumUnits() + unitIndex] = valueIn;
//		return 0;
//	} else {
//		return 1;
//	}
//}
//
//int phSolver::GetValue(SCField field, int unitIndex, int varIndex,
//		double& valueOut) {
//	//printField(field);
//	double *dataArray = (double*) UMap_FieldToPtr_[field];
//	if (NULL != dataArray) {
//		valueOut = dataArray[varIndex * field.GetNumUnits() + unitIndex];
//		return 0;
//	} else {
//		return 1;
//	}
//}
//
//int phSolver::SetValue(SCField field, int unitIndex, int varIndex,
//		int valueIn) {
//	//printField(field);
//	int *dataArray = (int*) UMap_FieldToPtr_[field];
//	if (NULL != dataArray) {
//		dataArray[varIndex * field.GetNumUnits() + unitIndex] = valueIn;
//		return 0;
//	} else {
//		return 1;
//	}
//}
//
//int phSolver::GetValue(SCField field, int unitIndex, int varIndex,
//		int& valueOut) {
//	//printField(field);
//	int *dataArray = (int*) UMap_FieldToPtr_[field];
//	if (NULL != dataArray) {
//		valueOut = dataArray[varIndex * field.GetNumUnits() + unitIndex];
//		return 0;
//	} else {
//		return 1;
//	}
//}

//int phSolver::GetNumBlocks() {
//	return Vec_BlockField_Ptr_.size();
//}
//
//int phSolver::GetBlockSize(int BlockId){
//	return Vec_BlockField_Field_[BlockId].GetNumUnits();
//}
//
//int phSolver::GetValueBlock(int BlockId, int unitIndex, int varIndex, int &value){
//
//	if (NULL != Vec_BlockField_Ptr_[BlockId]) {
//		value = Vec_BlockField_Ptr_[BlockId][varIndex * Vec_BlockField_Field_[BlockId].GetNumUnits() + unitIndex];
//		return 0;
//	} else {
//		return 1;
//	}
//}


int phSolver::getNumProcs() {
	return workfc.numpe;
}

int phSolver::getNumDofHolders() {
	return conpar.nshg;
}

int phSolver::getNumTimeSteps() {
	return inpdat.nstep[0];
}

int phSolver::getTime() {
	return timdat.lstep;
}

void phSolver::setTime(int timeIn) {
	timdat.lstep = timeIn; // this is dangerous to use right now
}

int phSolver::getSize() {

	int numVarsS = 4; // no scalar dof
	int numVarsD = NSD;

	int numUnits = conpar.nshguniq;

	return (numVarsD + 2 * numVarsS) * numUnits;
}

bool phSolver::hasFinished() {
	return timdat.istep >= this->getNumTimeSteps();
}

phSolver::~phSolver() {
}


void phSolverMessagePrinter::printDebugMsg(char* debugMsg) {
	//int myrank;
	//MPI_Comm_rank (this->iNewComm_C, &myrank);
	//if ( 0 == myrank ) {
	fprintf(stdout, "PHSOLVER_DEBUG - %s\n", debugMsg);
	//}
}

void phSolverMessagePrinter::printWarningMsg(char* warningMsg) {
	//int myrank;
	//MPI_Comm_rank (this->iNewComm_C, &myrank);
	//if ( 0 == myrank ) {
	fprintf(stderr, "PHSOLVER_WARNING - %s\n", warningMsg);
	//}
}

void phSolverMessagePrinter::printErrorMsg(char* errorMsg) {
	//int myrank;
	//MPI_Comm_rank (this->iNewComm_C, &myrank);
	//if ( 0 == myrank ) {
	fprintf(stderr, "PHSOLVER_ERROR - %s\n", errorMsg);
	//}
}
