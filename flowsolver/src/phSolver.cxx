/****************************************************************************** 

 (c) 2009-2010 Scientific Computation Research Center,
 Rensselaer Polytechnic Institute. All rights reserved.

 The LICENSE-SCOREC file included with this distribution describes the terms
 of the SCOREC Non-Commercial License this program is distributed under.

 *******************************************************************************/

#include "phSolver.h"

namespace phSolverMessagePrinter {
void print_error_code(int ierr);
void printDebugMsg(char* debugMsg);
void printErrorMsg(char* errorMsg);
void printWarningMsg(char* warningMsg);
}

void printField(SCField &inField) {
	printf(
			"name = %s, \
          domain type = %s, \
          number of units = %d, \
          number of variables = %d, \
          value type = %s\n",
			inField.GetName(),
			SCField::getDomainTypeName(inField.GetDomainType()),
			inField.GetNumUnits(), inField.GetNumVars(),
			SCField::getValTypeName(inField.GetValType()));
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

	/* Populate these lists based on what was read from solver.inp.
	 ex) 'LevelSet->Output level set gradient' set to 'True' will
	 add the graphi field to the list of produced fields.
	 The required fields list has a minimium number of entries; there
	 may be some solver.inp inputs that effect the list...
	 */
	initProducedFieldsList_(); // needed for when fields are not read from files
	initRequiredFieldsList_();
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

extern "C" void phSolverUpdateBlockField(char* pc_name, int i_numUnits, int i_numVars, int* pi_dataArray) {

    SCField inField(pc_name, "Volume", i_numUnits, i_numVars, "Integer");

	phSolver::Instance()->UpdateBlockField(inField, pi_dataArray);
}

void phSolver::UpdateBlockField(SCField inField, int* pi_dataArray) {
	if (NULL != pi_dataArray) {
		Vec_BlockField_Field_.push_back( inField );
		Vec_BlockField_Ptr_.push_back( pi_dataArray );
	}
}

/* 
 find the field in the list of required fields and remove it from the list - this function is meant to be called from common/readnblk.f and other locations where data is read from geombc and restart files
 */
extern "C" void phSolverRemoveField(char* pc_name, char* pc_domainType,
		char* pc_valType) {
	SCField inField(pc_name, pc_domainType, 0, 0, pc_valType);

	if (1 == phSolver::Instance()->RemoveField(inField)) {
		char warningMsg[256];
		sprintf(warningMsg, "field %s was not removed", inField.GetName());
		phSolverMessagePrinter::printWarningMsg(warningMsg);
	}
}

int phSolver::RemoveField(SCField inField) {
	bool i_notFound = 1;

	map<string,SCField>::iterator iter;


	for (iter = Vec_RequiredFields_.begin(); iter != Vec_RequiredFields_.end(); iter++) {
		if (0 == strcmp((*iter).second.GetName(), inField.GetName())
				&& (*iter).second.GetDomainType()
				== inField.GetDomainType()
				&& (*iter).second.GetValType()
				== inField.GetValType()) {
			i_notFound = 0;
			Vec_RequiredFields_.erase(iter);
			break;
		}
	}
	return i_notFound;
}

/* 
 find the field and update its data - this function is meant to be called from common/readnblk.f and other locations where data is read from geombc and restart files
 */
// see if it is possible to avoid using a global function here
extern "C" void phSolverUpdateField(char* pc_name, char* pc_domainType,
		char* pc_valType, int i_numUnits, int i_numVar, double *pd_dataArray,
		int *pi_dataArray) {
	SCField inField(pc_name, pc_domainType, i_numUnits, i_numVar, pc_valType);

//	char debugMessage[256];
//	sprintf(debugMessage, "pointer to field %s: dbl: %d int: %i",
//			inField.GetName(), pd_dataArray, pi_dataArray);
//	phSolverMessagePrinter::printDebugMsg(debugMessage);

	void* pv_dataArray = NULL;
	if (pd_dataArray != NULL && pi_dataArray == NULL) {
		pv_dataArray = (void*) pd_dataArray;
	} else if (pd_dataArray == NULL && pi_dataArray != NULL) {
		pv_dataArray = (void*) pi_dataArray;
	} else {
//		char warningMsg[256];
//		sprintf(warningMsg, "pointer to field %s was not updated",
//				inField.GetName());
//		phSolverMessagePrinter::printWarningMsg(warningMsg);
	}

	if (1 == phSolver::Instance()->UpdateField(inField, pv_dataArray)) {
//		char warningMsg[256];
//		sprintf(warningMsg, "field %s was not updated", inField.GetName());
//		phSolverMessagePrinter::printWarningMsg(warningMsg);
	}

}

int phSolver::UpdateField(SCField inField, void* pv_dataArray) {
	bool i_notFound = 1;

	map<string,SCField>::iterator iter;

	for (iter = Vec_RequiredFields_.begin(); iter != Vec_RequiredFields_.end(); iter++) {
		if (0 == strcmp((*iter).second.GetName(), inField.GetName())
				&& (*iter).second.GetDomainType()
				== inField.GetDomainType()
				&& (*iter).second.GetValType()
				== inField.GetValType()) {
			(*iter).second.SetNumUnits(inField.GetNumUnits());
			(*iter).second.SetNumVars(inField.GetNumVars());
			i_notFound = 0;
			break;
		}
	}

	if (0 == i_notFound) {
		Vec_SatisfiedFields_.insert( pair<string,SCField>(string(inField.GetName()),inField) );
		if (NULL != pv_dataArray) {
			UMap_FieldToPtr_.insert( pair<SCField,void*>(inField,pv_dataArray) );
		}
	}

	return i_notFound;
}

// when to set these pointers??
int phSolver::initProducedFieldsList_() {

	return 0;
}

int phSolver::NumRequiredFields() const {
	return Vec_RequiredFields_.size();
}

const SCField* phSolver::GetRequiredField(string fieldIdx) {

	return &(Vec_RequiredFields_[fieldIdx]);

}

//const SCField* phSolver::GetProducedField(string fieldIdx) {
//
//	return &(Vec_ProducedFields_[fieldIdx]);
//
//}

//const SCField phSolver::GetField(SCField infield, int* err) {
//	SCField matchingField;
//	*err = 1;
//	for (int i = 0; i < Vec_ProducedFields_.size(); i++) {
//		if (infield == Vec_ProducedFields_[i]) {
//			matchingField = Vec_ProducedFields_[i];
//			*err = 0;
//			break;
//		}
//	}
//	return matchingField;
//}

int phSolver::initRequiredFieldsList_() {

	MPI_Comm_size(this->iNewComm_C, &numProcs_);

	//HACK hardcoded
	//HACK set numDofs to '-1' to indicate that it is not initialized - until it is initialized it prevents the class comparison operator from being used

	SCField temp;

	// from restart.*.*
	Vec_RequiredFields_.insert( pair<string,SCField>("solution",
			SCField("solution", SCField::Volume, -1, 7, SCField::Double)) );

	Vec_RequiredFields_.insert( pair<string,SCField>("time derivative of solution",
			SCField("time derivative of solution", SCField::Volume, -1, 7, SCField::Double)) ); // rnb

	if (nomodule.ideformwall > 0) {
		Vec_RequiredFields_.insert( pair<string,SCField>("displacement",
				SCField("displacement", SCField::Volume, -1, 3, SCField::Double)) ); //rnb - this is only required if CARDIOVASCULAR parameters in solver.inp is enabled
		Vec_RequiredFields_.insert( pair<string,SCField>("elastic modulus scalar",
				SCField("elastic modulus scalar", SCField::Scalar, 1, 1, SCField::Double)) ); //vlm - this is only required if CARDIOVASCULAR parameters in solver.inp is enabled

		Vec_RequiredFields_.insert( pair<string,SCField>("observation function displacement",
					SCField("observation function displacement", SCField::Volume, -1, 7, SCField::Integer)) );
	}

	Vec_RequiredFields_.insert( pair<string,SCField>("observation function solution",
			SCField("observation function solution", SCField::Volume, -1, 7, SCField::Integer)) );

	Vec_RequiredFields_.insert( pair<string,SCField>("observation function time derivative of solution",
			SCField("observation function time derivative of solution", SCField::Volume, -1, 7, SCField::Integer)) );

	Vec_RequiredFields_.insert( pair<string,SCField>("temporary_array",
				SCField("temporary_array", SCField::Volume, -1, 7, SCField::Double)) ); //

	// from geombc.dat.*
	Vec_RequiredFields_.insert( pair<string,SCField>("number of nodes",
			SCField("number of nodes", SCField::Scalar, 1, 1,SCField::Integer)) ); // rnb
	Vec_RequiredFields_.insert( pair<string,SCField>("number of modes",
			SCField("number of modes", SCField::Scalar, 1, 1,SCField::Integer)) ); // rnb
	Vec_RequiredFields_.insert( pair<string,SCField>("number of shape functions solved on processor",
			SCField("number of shape functions solved on processor",SCField::Scalar, 1, 1, SCField::Integer)) ); // the geombc.dat file this is based on has typos in the string // rnb
	Vec_RequiredFields_.insert( pair<string,SCField>("number of interior elements",
			SCField("number of interior elements", SCField::Scalar, 1, 1,SCField::Integer)) ); // rnb
	Vec_RequiredFields_.insert( pair<string,SCField>("number of boundary elements",
			SCField("number of boundary elements", SCField::Scalar, 1, 1,SCField::Integer)) ); // rnb
	Vec_RequiredFields_.insert( pair<string,SCField>("maximum number of element nodes",
			SCField("maximum number of element nodes", SCField::Scalar, 1, 1, SCField::Integer)) ); // rnb
	Vec_RequiredFields_.insert( pair<string,SCField>("number of interior tpblocks",
			SCField("number of interior tpblocks", SCField::Scalar, 1, 1,SCField::Integer)) ); // rnb - number of types of blocks on interior
	Vec_RequiredFields_.insert( pair<string,SCField>("number of boundary tpblocks",
			SCField("number of boundary tpblocks", SCField::Scalar, 1, 1,SCField::Integer)) ); // rnb -  number of types of blocks on boundary
	Vec_RequiredFields_.insert( pair<string,SCField>("number of nodes with Dirichlet BCs",
			SCField("number of nodes with Dirichlet BCs", SCField::Scalar, 1, 1,SCField::Integer)) ); // rnb
	Vec_RequiredFields_.insert( pair<string,SCField>("co-ordinates",
			SCField("co-ordinates", SCField::Volume, -1, 3, SCField::Double)) ); // rnb
//		Vec_RequiredFields_.push_back(
//				SCField("connectivity interior linear tetrahedron", SCField::Volume,
//						-1, -1, SCField::Integer)); // genblk
//		Vec_RequiredFields_.push_back(
//				SCField("connectivity boundary linear tetrahedron",
//						SCField::Surface, -1, -1, SCField::Integer)); // genblkb
//		Vec_RequiredFields_.push_back(
//				SCField("nbc codes linear tetrahedron", SCField::Surface, -1, -1,
//						SCField::Integer)); // genblkb
//		Vec_RequiredFields_.push_back(
//				SCField("nbc values linear tetrahedron", SCField::Surface, -1, -1,
//						SCField::Double)); // genblkb
	Vec_RequiredFields_.insert( pair<string,SCField>("bc mapping array",
			SCField("bc mapping array", SCField::Volume, -1, 1,SCField::Integer)) ); // rnb
	Vec_RequiredFields_.insert( pair<string,SCField>("bc codes array",
			SCField("bc codes array", SCField::Volume, -1, 1,SCField::Integer)) ); // rnb
	Vec_RequiredFields_.insert( pair<string,SCField>("boundary condition array",
			SCField("boundary condition array", SCField::Volume, -1, 1,SCField::Integer)) ); // rnb
	//Vec_RequiredFields_.push_back(SCField("periodic masters array", SCField::Volume, -1, 1, SCField::Integer));  // rnb
	// add the communication fields if running in parallel
	if (numProcs_ > 1) {
		Vec_RequiredFields_.insert( pair<string,SCField>("size of ilwork array",
				SCField("size of ilwork array", SCField::Scalar, 1, 1,SCField::Integer)) ); // rnb
		Vec_RequiredFields_.insert( pair<string,SCField>("ilwork",
				SCField("ilwork", SCField::Volume, -1, 1, SCField::Integer)) ); // rnb
	}
	//  Vec_RequiredFields_.push_back(SCField("number of nodes in the mesh", SCField::Scalar, 1, 1, SCField::Integer)); // NO - read in partition.cc for mesh partitioning - not required
	//  Vec_RequiredFields_.push_back(SCField("number of edges in the mesh", SCField::Scalar, 1, 1, SCField::Integer)); // NO - read in partition.cc for mesh partitioning - not required
	//  Vec_RequiredFields_.push_back(SCField("number of faces in the mesh", SCField::Scalar, 1, 1, SCField::Integer)); // NO - not read from file upon restart
	//  Vec_RequiredFields_.push_back(SCField("mode number map from partition to global", SCField::Volume, -1, 1, SCField::Integer)); // NO  - read in partition.cc for mesh partitioning - not required
	//  Vec_RequiredFields_.push_back(SCField("number of processors", SCField::Scalar, 1, 1, SCField::Integer)); // NO - read in partition.cc for mesh partitioning - not required
	//  Vec_RequiredFields_.push_back(SCField("number of global modes", SCField::Scalar, 1, 1, SCField::Integer));  // NO - read in partition.cc for mesh partitioning - not required

	//cout << "total of " << Vec_RequiredFields_.size() << " required fields" << endl;

	Vec_RequiredFields_.insert( pair<string,SCField>("local index of unique nodes",
			SCField("local index of unique nodes", SCField::Volume, -1, 1,SCField::Integer)) ); // rnb

	return 0;
}

int phSolver::SolverInit() {
	int myrank;
	MPI_Comm_rank(this->iNewComm_C, &myrank);

	char errorMsg[256], warningMsg[256];

	int i_numReqFieldsNotSatisfied = 0;
	map<string,SCField>::iterator fieldIter;

	if (-1 == chdir(pc_pathToProcsCaseDir_)) {
		sprintf(errorMsg, "chdir(%s) failed", pc_pathToProcsCaseDir_);
		phSolverMessagePrinter::printErrorMsg(errorMsg);
		return 1;
	}
//	phSolverMessagePrinter::printDebugMsg("Initializing solver...");
	itrdrv_init(); // calls the fortran routine

//	for (fieldIter = Vec_RequiredFields_.begin();
//			fieldIter != Vec_RequiredFields_.end(); fieldIter++) {
//		map<string,SCField>::iterator satFieldIter;
//		satFieldIter = find(Vec_SatisfiedFields_.begin(),
//				Vec_SatisfiedFields_.end(), *fieldIter);
//		if (satFieldIter == Vec_SatisfiedFields_.end()) {
//			sprintf(warningMsg, "field %s is missing", fieldIter->GetName());
//			phSolverMessagePrinter::printWarningMsg(warningMsg);
//			++i_numReqFieldsNotSatisfied;
//		}
//	}
//	if (0 == i_numReqFieldsNotSatisfied) {
//
//		return 0;
//	} else {
//		sprintf(errorMsg,
//				"solver cannot run because all required fields have not been satisfied");
//		phSolverMessagePrinter::printErrorMsg(errorMsg);
//		return 1;
//	}

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

int phSolver::SetValue(SCField field, int unitIndex, int varIndex,
		double valueIn) {
	//printField(field);
	double *dataArray = (double*) UMap_FieldToPtr_[field];
	if (NULL != dataArray) {
		dataArray[varIndex * field.GetNumUnits() + unitIndex] = valueIn;
		return 0;
	} else {
		return 1;
	}
}

int phSolver::GetValue(SCField field, int unitIndex, int varIndex,
		double& valueOut) {
	//printField(field);
	double *dataArray = (double*) UMap_FieldToPtr_[field];
	if (NULL != dataArray) {
		valueOut = dataArray[varIndex * field.GetNumUnits() + unitIndex];
		return 0;
	} else {
		return 1;
	}
}

int phSolver::SetValue(SCField field, int unitIndex, int varIndex,
		int valueIn) {
	//printField(field);
	int *dataArray = (int*) UMap_FieldToPtr_[field];
	if (NULL != dataArray) {
		dataArray[varIndex * field.GetNumUnits() + unitIndex] = valueIn;
		return 0;
	} else {
		return 1;
	}
}

int phSolver::GetValue(SCField field, int unitIndex, int varIndex,
		int& valueOut) {
	//printField(field);
	int *dataArray = (int*) UMap_FieldToPtr_[field];
	if (NULL != dataArray) {
		valueOut = dataArray[varIndex * field.GetNumUnits() + unitIndex];
		return 0;
	} else {
		return 1;
	}
}

int phSolver::GetNumBlocks() {
	return Vec_BlockField_Ptr_.size();
}

int phSolver::GetBlockSize(int BlockId){
	return Vec_BlockField_Field_[BlockId].GetNumUnits();
}

int phSolver::GetValueBlock(int BlockId, int unitIndex, int varIndex, int &value){

	if (NULL != Vec_BlockField_Ptr_[BlockId]) {
		value = Vec_BlockField_Ptr_[BlockId][varIndex * Vec_BlockField_Field_[BlockId].GetNumUnits() + unitIndex];
		return 0;
	} else {
		return 1;
	}
}


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

	int numVarsS = this->GetRequiredField("solution")->GetNumVars() - 1; // get rid of scalar dof
	int numVarsD = this->GetRequiredField("displacement")->GetNumVars();

	int numUnits = this->GetRequiredField("local index of unique nodes")->GetNumUnits();

	return (numVarsD + 2 * numVarsS) * numUnits;
}

bool phSolver::hasFinished() {
	return timdat.istep >= this->getNumTimeSteps();
}

int phSolver::input_fform() {

	int ierr = 0;
	int i, j;
	char* path_to_config = 0, *inpfilename_env = 0;
	char complete_filename[256], inpfname[256];

	try {
		// get the input file stream
		path_to_config = getenv("PHASTA_CONFIG");
		if (path_to_config)
			strcpy(complete_filename, path_to_config);
		else
			strcpy(complete_filename, ".");
		strcat(complete_filename, "/input.config");
		inpfilename_env = getenv("PHASTA_INPFILE");
		if (inpfilename_env)
			strcpy(inpfname, inpfilename_env);
		else
			strcpy(inpfname, "solver.inp");
		if (workfc.myrank == workfc.master) {
			printf("\n Complete Filename: %s \n", complete_filename);
			printf("\n Local Config: %s \n\n", inpfname);
		}
		printf("\n Complete Filename: %s \n", complete_filename);
		printf("\n Local Config: %s \n\n", inpfname);
		string def(complete_filename);
		CInput inp(inpfname, def);

		// Disabled Features

		conpar.iALE = inp.GetValue("iALE");
		conpar.icoord = inp.GetValue("icoord");
		conpar.irs = inp.GetValue("irs");
		conpar.iexec = inp.GetValue("iexec");
		timpar.ntseq = inp.GetValue("ntseq");
		solpar.imap = inp.GetValue("imap");

		// Solution Control Keywords

		if ((string) inp.GetValue("Equation of State") == "Incompressible")
			matdat.matflg[0][0] = -1;
		if ((string) inp.GetValue("Equation of State") == "Compressible")
			matdat.matflg[0][0] = 0;
		inpdat.Delt[0] = inp.GetValue("Time Step Size");
		inpdat.nstep[0] = inp.GetValue("Number of Timesteps");
		if ((string) inp.GetValue("Viscous Control") == "Viscous")
			conpar.navier = 1;
		else
			conpar.navier = 0;

		if ((string) inp.GetValue("Turbulence Model") == "No-Model") {
			turbvari.irans = 0;
			turbvari.iles = 0;
		} else if ((string) inp.GetValue("Turbulence Model") == "LES") {
			turbvari.iles = 1;
			turbvari.irans = 0;
		} else if ((string) inp.GetValue("Turbulence Model") == "RANS-SA") {
			turbvari.iles = 0;
			turbvari.irans = -1;
		} else if ((string) inp.GetValue("Turbulence Model") == "RANS") {
			turbvari.iles = 0;
			turbvari.irans = -1; // assume S-A for backward compatibility
		} else if ((string) inp.GetValue("Turbulence Model") == "RANS-KE") {
			turbvari.iles = 0;
			turbvari.irans = -2;
		} else if ((string) inp.GetValue("Turbulence Model") == "DES") {
			turbvari.iles = 1;
			turbvari.irans = -1;
		} else {
			cout
					<< " Turbulence Model: Only Legal Values ( No-Model, LES, RANS-SA, RANS-KE, DES )";
			cout << endl;
			exit(1);
		}

		if (turbvari.iles * turbvari.irans != 0)
			turbvar.eles = inp.GetValue("DES Edge Length");

		int solflow, solheat, solscalr, ilset;
		((string) inp.GetValue("Solve Flow") == "True") ?
				solflow = 1 : solflow = 0;
		((string) inp.GetValue("Solve Heat") == "True") ?
				solheat = 1 : solheat = 0;
		//for compressible solheat= False so
		if ((string) inp.GetValue("Equation of State") == "Compressible")
			solheat = 0;
		ilset = (int) inp.GetValue("Solve Level Set");
		solscalr = (int) inp.GetValue("Solve Scalars");
		solscalr += ilset;
		if (turbvari.irans == -1)
			solscalr++;
		if (turbvari.irans == -2)
			solscalr = solscalr + 2;
		if (solscalr > 4) {
			cout << " Only Four Scalars are supported \n";
			cout << " Please reduce number of scalars \n";
			exit(1);
		}
		inpdat.impl[0] = 10 * solflow + solscalr * 100 + solheat;

		levlset.iLSet = ilset;
		if (ilset > 0) {
			levlset.epsilon_ls = inp.GetValue(
					"Number of Elements Across Interface");
			levlset.epsilon_lsd = inp.GetValue(
					"Number of Elements Across Interface for Redistancing");
			levlset.dtlset = inp.GetValue("Pseudo Time step for Redistancing");
			levlset.iExpLSSclr2 = inp.GetValue(
					"Explicit Solve for Redistance Field");
			levlset.iExpLSSclr1 = inp.GetValue(
					"Explicit Solve for Scalar 1 Field");
			if ((string) inp.GetValue("Apply Volume Constraint") == "True") {
				levlset.ivconstraint = 1;
			} else if ((string) inp.GetValue("Apply Volume Constraint")
					== "False") {
				levlset.ivconstraint = 0;
			} else {
				cout
						<< "Apply Volume Constraint: Only Legal Values (True, False) ";
				cout << endl;
				exit(1);
			}
		}

		vector<double> vec;

		// OUTPUT CONTROL KEY WORDS.

		conpar.necho = inp.GetValue("Verbosity Level");
		outpar.ntout = inp.GetValue("Number of Timesteps between Restarts");
		if ((string) inp.GetValue("Print Statistics") == "True")
			outpar.ioform = 2;
		else
			outpar.ioform = 1;

		if ((string) inp.GetValue("Print Wall Fluxes") == "True")
			outpar.iowflux = 1;
		else
			outpar.iowflux = 0;

		if ((string) inp.GetValue("Print FieldView") == "True")
			outpar.iofieldv = 1;
		else
			outpar.iofieldv = 0;

		if ((string) inp.GetValue("Print ybar") == "True")
			outpar.ioybar = 1;
		else
			outpar.ioybar = 0;

		strcpy(outpar.iotype,
				((string) inp.GetValue("Data Block Format")).c_str());
		//strcpy( phasta_iotype , ((string)inp.GetValue("Data Block Format")).c_str());
		turbvar.sonfathvar = inp.GetValue("Number of Father Nodes");

		if ((string) inp.GetValue("Print Residual at End of Step") == "True")
			genpar.lstres = 1;
		else
			genpar.lstres = 0;

		if ((string) inp.GetValue("Print Error Indicators") == "True")
			turbvar.ierrcalc = 1;
		else
			turbvar.ierrcalc = 0;

		if ((string) inp.GetValue("Print Velocity Hessian") == "True")
			turbvar.ihessian = 1;
		else
			turbvar.ihessian = 0;

		if (turbvar.ierrcalc == 1)
			turbvari.ierrsmooth = inp.GetValue(
					"Number of Error Smoothing Iterations");

		int nsrfCM = inp.GetValue("Number of Force Surfaces");
		if (nsrfCM > 0) {
			vector<int> ivec = inp.GetValue(
					"Surface ID's for Force Calculation");
			for (i = 0; i < MAXSURF + 1; i++)
				aerfrc.nsrflist[i] = 0;
			for (i = 0; i < nsrfCM; i++) {
				aerfrc.nsrflist[ivec[i]] = 1;
				//        cout <<"surface in force list "<< ivec[i] << endl;
			}
			ivec.erase(ivec.begin(), ivec.end());
		}

		aerfrc.isrfIM = inp.GetValue("Surface ID for Integrated Mass");
		//Limiting
		vec = inp.GetValue("Limit u1");
		for (i = 0; i < 3; i++) {
			turbvar.ylimit[0][i] = vec[i];
		}
		vec.erase(vec.begin(), vec.end());

		vec = inp.GetValue("Limit u2");
		for (i = 0; i < 3; i++) {
			turbvar.ylimit[1][i] = vec[i];
		}
		vec.erase(vec.begin(), vec.end());

		vec = inp.GetValue("Limit u3");
		for (i = 0; i < 3; i++) {
			turbvar.ylimit[2][i] = vec[i];
		}
		vec.erase(vec.begin(), vec.end());

		vec = inp.GetValue("Limit Pressure");
		for (i = 0; i < 3; i++) {
			turbvar.ylimit[3][i] = vec[i];
		}
		vec.erase(vec.begin(), vec.end());

		vec = inp.GetValue("Limit Temperature");
		for (i = 0; i < 3; i++) {
			turbvar.ylimit[4][i] = vec[i];
		}
		vec.erase(vec.begin(), vec.end());

		//Material Properties Keywords
		matdat.nummat = levlset.iLSet + 1;
		if ((string) inp.GetValue("Shear Law") == "Constant Viscosity")
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][1] = 0;

		if ((string) inp.GetValue("Bulk Viscosity Law")
				== "Constant Bulk Viscosity")
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][2] = 0;

		mmatpar.pr = inp.GetValue("Prandtl Number");

		if ((string) inp.GetValue("Conductivity Law")
				== "Constant Conductivity")
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][3] = 0;

		vec = inp.GetValue("Density");
		for (i = 0; i < levlset.iLSet + 1; i++) {
			matdat.datmat[i][0][0] = vec[i];
		}
		vec.erase(vec.begin(), vec.end());

		vec = inp.GetValue("Viscosity");
		for (i = 0; i < levlset.iLSet + 1; i++) {
			matdat.datmat[i][1][0] = vec[i];
		}
		vec.erase(vec.begin(), vec.end());

		//      vec = inp.GetValue("Specific Heat");
		for (i = 0; i < levlset.iLSet + 1; i++) {
			matdat.datmat[i][2][0] = 0;
		}
		//      vec.erase(vec.begin(),vec.end());

		vec = inp.GetValue("Thermal Conductivity");
		for (i = 0; i < levlset.iLSet + 1; i++) {
			matdat.datmat[i][3][0] = vec[i];
		}
		vec.erase(vec.begin(), vec.end());

		vec = inp.GetValue("Scalar Diffusivity");
		for (i = 0; i < solscalr; i++) {
			sclrs.scdiff[i] = vec[i];
		}
		vec.erase(vec.begin(), vec.end());

		if ((string) inp.GetValue("Zero Mean Pressure") == "True")
			turbvar.pzero = 1;

		turbvar.rmutarget = inp.GetValue("Target Viscosity For Step NSTEP");

		if ((string) inp.GetValue("Body Force Option") == "None") {
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][4] = 0;
		} else if ((string) inp.GetValue("Body Force Option") == "Vector") {
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][4] = 1;
		} else if ((string) inp.GetValue("Body Force Option")
				== "User e3source.f") {
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][4] = 3;
		} else if ((string) inp.GetValue("Body Force Option") == "Boussinesq") {
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][4] = 2;
		} else if ((string) inp.GetValue("Body Force Option")
				== "Cooling Analytic") {
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][4] = 4;
		} else if ((string) inp.GetValue("Body Force Option")
				== "Cooling Initial Condition") {
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][4] = 5;
		}

		// the following block of stuff is common to all cooling type sponges.
		// Specific stuff belongs in the conditionals above

		if (matdat.matflg[0][4] >= 4) {
			spongevar.betamax = inp.GetValue(
					"Maximum Value of Sponge Parameter");
			spongevar.zinsponge = inp.GetValue(
					"Inflow Cooling Sponge Ends at z");
			spongevar.zoutsponge = inp.GetValue(
					"Outflow Cooling Sponge Begins at z");
			spongevar.radsponge = inp.GetValue(
					"Radial Cooling Sponge Begins at r");
			spongevar.grthosponge = inp.GetValue(
					"Sponge Growth Coefficient Outflow");
			spongevar.grthisponge = inp.GetValue(
					"Sponge Growth Coefficient Inflow");

			spongevar.spongecontinuity = 0;
			spongevar.spongemomentum1 = 0;
			spongevar.spongemomentum2 = 0;
			spongevar.spongemomentum3 = 0;
			spongevar.spongeenergy = 0;

			if ((string) inp.GetValue("Sponge for Continuity Equation")
					== "True")
				spongevar.spongecontinuity = 1;
			if ((string) inp.GetValue("Sponge for x Momentum Equation")
					== "True")
				spongevar.spongemomentum1 = 1;
			if ((string) inp.GetValue("Sponge for y Momentum Equation")
					== "True")
				spongevar.spongemomentum2 = 1;
			if ((string) inp.GetValue("Sponge for z Momentum Equation")
					== "True")
				spongevar.spongemomentum3 = 1;
			if ((string) inp.GetValue("Sponge for Energy Equation") == "True")
				spongevar.spongeenergy = 1;

		}

		vec = inp.GetValue("Body Force");
		for (i = 0; i < levlset.iLSet + 1; i++) {
			matdat.datmat[i][4][0] = vec[0 + i * 3];
			matdat.datmat[i][4][1] = vec[1 + i * 3];
			matdat.datmat[i][4][2] = vec[2 + i * 3];
		}
		vec.erase(vec.begin(), vec.end());

		vec = inp.GetValue("Body Force Pressure Gradient");
		for (i = 0; i < levlset.iLSet + 1; i++) {
			matdat.datmat[i][6][0] = vec[0 + i * 3];
			matdat.datmat[i][6][1] = vec[1 + i * 3];
			matdat.datmat[i][6][2] = vec[2 + i * 3];
		}
		vec.erase(vec.begin(), vec.end());

		if ((string) inp.GetValue("Surface Tension Option") == "No") {
			genpar.isurf = 0;
		} else if ((string) inp.GetValue("Surface Tension Option") == "Yes") {
			genpar.isurf = 1;
		} else {
			cout << " Surface Tension: Only Legal Values (Yes, No) ";
			cout << endl;
			exit(1);
		}
		if (genpar.isurf > 0) {
			genpar.Bo = inp.GetValue("Bond Number");
		}

		genpar.EntropyPressure = inp.GetValue(
				"Entropy Form of Pressure Constraint on Weight Space");

		if ((string) inp.GetValue("Rotating Frame of Reference") == "True") {
			matdat.matflg[0][5] = 1;
			vec = inp.GetValue("Rotating Frame of Reference Rotation Rate");
			matdat.datmat[0][5][0] = vec[0];
			matdat.datmat[0][5][1] = vec[1];
			matdat.datmat[0][5][2] = vec[2];
			vec.erase(vec.begin(), vec.end());
		} else {
			matdat.matflg[0][5] = 0;
			matdat.datmat[0][5][0] = 0.;
			matdat.datmat[0][5][1] = 0.;
			matdat.datmat[0][5][2] = 0.;
		}

		//Linear Solver parameters
		if ((string) inp.GetValue("Solver Type")
				== "ACUSIM with P Projection") {
			incomp.iprjFlag = 0;
			incomp.ipresPrjFlag = 1;
		} else if ((string) inp.GetValue("Solver Type") == "ACUSIM") {
			incomp.iprjFlag = 0;
			incomp.ipresPrjFlag = 0;
		} else if ((string) inp.GetValue("Solver Type")
				== "ACUSIM with Velocity Projection") {
			incomp.iprjFlag = 1;
			incomp.ipresPrjFlag = 0;
		} else if ((string) inp.GetValue("Solver Type")
				== "ACUSIM with Full Projection") {
			incomp.iprjFlag = 1;
			incomp.ipresPrjFlag = 1;
		} else if ((string) inp.GetValue("Solver Type")
				== "GMRES Matrix Free") {
			inpdat.impl[0] += 10 * solflow;
		} else if ((string) inp.GetValue("Solver Type") == "GMRES EBE") {
			inpdat.impl[0] += 20 * solflow;
		}
		//GMRES sparse is assumed default and has the value of 10, MFG 20,
		// EBE 30

		//    inpdat.niter[0] = inp.GetValue("Number of Solves per Time Step");
		solpar.nGMRES = inp.GetValue("Number of GMRES Sweeps per Solve");
		solpar.Kspace = inp.GetValue(
				"Number of Krylov Vectors per GMRES Sweep");
		inpdat.LHSupd[0] = inp.GetValue(
				"Number of Solves per Left-hand-side Formation");
		inpdat.epstol[0] = inp.GetValue("Tolerance on Momentum Equations");
		incomp.prestol = inp.GetValue(
				"Tolerance on ACUSIM Pressure Projection");
		incomp.minIters = inp.GetValue(
				"Minimum Number of Iterations per Nonlinear Iteration");
		incomp.maxIters = inp.GetValue(
				"Maximum Number of Iterations per Nonlinear Iteration");
		inpdat.deltol[0][0] = inp.GetValue("Velocity Delta Ratio");
		inpdat.deltol[1][0] = inp.GetValue("Pressure Delta Ratio");
		incomp.nPrjs = inp.GetValue("Number of Velocity Projection Vectors");
		incomp.nPresPrjs = inp.GetValue(
				"Number of Pressure Projection Vectors");
		incomp.iverbose = inp.GetValue("ACUSIM Verbosity Level");

		if (solheat == 1) {
			inpdat.epstol[1] = inp.GetValue("Temperature Solver Tolerance");
			inpdat.LHSupd[1] =
					inp.GetValue(
							"Number of Solves of Temperature per Left-hand-side Formation");
		}

		// The following is where you should put any inputs that are able to
		// input differently for each scalar.  It is a little tedious in the code
		// but it should make the solver.inp easier to understand. Note this will
		// require some care with regression tests.

		if (solscalr > 0) {
			inpdat.epstol[2] = inp.GetValue("Scalar 1 Solver Tolerance");
			inpdat.LHSupd[2] =
					inp.GetValue(
							"Number of Solves of Scalar 1 per Left-hand-side Formation");

			vec = inp.GetValue("Limit Scalar 1");
			for (i = 0; i < 3; i++) {
				turbvar.ylimit[5][i] = vec[i];
			}
			vec.erase(vec.begin(), vec.end());
		}

		if (solscalr > 1) {
			inpdat.epstol[3] = inp.GetValue("Scalar 2 Solver Tolerance");
			inpdat.LHSupd[3] =
					inp.GetValue(
							"Number of Solves of Scalar 2 per Left-hand-side Formation");

			vec = inp.GetValue("Limit Scalar 2");
			for (i = 0; i < 3; i++) {
				turbvar.ylimit[6][i] = vec[i];
			}
			vec.erase(vec.begin(), vec.end());
		}

		if (solscalr > 2) {
			inpdat.epstol[4] = inp.GetValue("Scalar 3 Solver Tolerance");
			inpdat.LHSupd[4] =
					inp.GetValue(
							"Number of Solves of Scalar 3 per Left-hand-side Formation");

			vec = inp.GetValue("Limit Scalar 3");
			for (i = 0; i < 3; i++) {
				turbvar.ylimit[7][i] = vec[i];
			}
			vec.erase(vec.begin(), vec.end());
		}

		if (solscalr > 3) {
			inpdat.epstol[5] = inp.GetValue("Scalar 4 Solver Tolerance");
			inpdat.LHSupd[5] =
					inp.GetValue(
							"Number of Solves of Scalar 4 per Left-hand-side Formation");

			vec = inp.GetValue("Limit Scalar 4");
			for (i = 0; i < 3; i++) {
				turbvar.ylimit[8][i] = vec[i];
			}
			vec.erase(vec.begin(), vec.end());
		}

		// DISCRETIZATION CONTROL

		genpar.ipord = inp.GetValue("Basis Function Order");
		if ((string) inp.GetValue("Time Integration Rule") == "First Order")
			inpdat.rhoinf[0] = -1;
		else
			inpdat.rhoinf[0] = (double) inp.GetValue(
					"Time Integration Rho Infinity");
		if ((string) inp.GetValue("Predictor at Start of Step")
				== "Same Velocity")
			genpar.ipred = 1;
		if ((string) inp.GetValue("Predictor at Start of Step")
				== "Zero Acceleration")
			genpar.ipred = 2;
		if ((string) inp.GetValue("Predictor at Start of Step")
				== "Same Acceleration")
			genpar.ipred = 3;
		if ((string) inp.GetValue("Predictor at Start of Step") == "Same Delta")
			genpar.ipred = 4;

		if ((string) inp.GetValue("Weak Form") == "Galerkin")
			solpar.ivart = 1;
		if ((string) inp.GetValue("Weak Form") == "SUPG")
			solpar.ivart = 2;

		if ((string) inp.GetValue("Flow Advection Form") == "Convective")
			solpar.iconvflow = 2;
		else if ((string) inp.GetValue("Flow Advection Form") == "Conservative")
			solpar.iconvflow = 1;
		if ((string) inp.GetValue("Scalar Advection Form") == "Convective")
			solpar.iconvsclr = 2;
		else if ((string) inp.GetValue("Scalar Advection Form")
				== "Conservative")
			solpar.iconvsclr = 1;
		if ((string) inp.GetValue("Use Conservative Scalar Convection Velocity")
				== "True")
			sclrs.consrv_sclr_conv_vel = 1;
		else if ((string) inp.GetValue(
				"Use Conservative Scalar Convection Velocity") == "False")
			sclrs.consrv_sclr_conv_vel = 0;
		// TAU INPUT
		if ((string) inp.GetValue("Tau Matrix") == "Diagonal-Shakib")
			genpar.itau = 0;
		else if ((string) inp.GetValue("Tau Matrix") == "Diagonal-Franca")
			genpar.itau = 1;
		else if ((string) inp.GetValue("Tau Matrix") == "Diagonal-Jansen(dev)")
			genpar.itau = 2;
		else if ((string) inp.GetValue("Tau Matrix") == "Diagonal-Compressible")
			genpar.itau = 3;
		else if ((string) inp.GetValue("Tau Matrix") == "Matrix-Mallet")
			genpar.itau = 10;
		else if ((string) inp.GetValue("Tau Matrix") == "Matrix-Modal")
			genpar.itau = 11;

		genpar.dtsfct = inp.GetValue("Tau Time Constant");
		genpar.taucfct = inp.GetValue("Tau C Scale Factor");

		// FLOW DISCONTINUITY CAPTURING

		if ((string) inp.GetValue("Discontinuity Capturing") == "Off")
			solpar.iDC = 0;
		else if ((string) inp.GetValue("Discontinuity Capturing")
				== "DC-mallet")
			solpar.iDC = 1;
		else if ((string) inp.GetValue("Discontinuity Capturing")
				== "DC-quadratic")
			solpar.iDC = 2;
		else if ((string) inp.GetValue("Discontinuity Capturing")
				== "DC-minimum")
			solpar.iDC = 3;
		else {
			cout << "Condition not defined for Discontinuity Capturing \n ";
			exit(1);
		}

		// SCALAR DISCONTINUITY CAPTURING

		vector<int> ivec = inp.GetValue("Scalar Discontinuity Capturing");
		for (i = 0; i < 2; i++)
			solpar.idcsclr[i] = ivec[i];
		ivec.erase(ivec.begin(), ivec.end());

		//        if((string)inp.GetValue("Scalar Discontinuity Capturing") == "No") solpar.idcsclr = 0;
		//      else if((string)inp.GetValue("Scalar Discontinuity Capturing") == "1") solpar.idcsclr = 1;
		//   else if((string)inp.GetValue("Scalar Discontinuity Capturing") == "2") solpar.idcsclr = 2;
		//   else {
		//        cout<< "Condition not defined for Scalar Discontinuity Capturing \n ";
		//        exit(1);
		//      }
		if ((string) inp.GetValue("Include Viscous Correction in Stabilization")
				== "True") {
			if (genpar.ipord == 1)
				genpar.idiff = 1;
			else
				genpar.idiff = 2;
		} else {
			genpar.idiff = 0;
		}

		timdat.flmpl = inp.GetValue("Lumped Mass Fraction on Left-hand-side");
		timdat.flmpr = inp.GetValue("Lumped Mass Fraction on Right-hand-side");

		if ((string) inp.GetValue("Dump CFL") == "True")
			timdat.iCFLworst = 1;

		intdat.intg[0][0] = inp.GetValue("Quadrature Rule on Interior");
		intdat.intg[0][1] = inp.GetValue("Quadrature Rule on Boundary");
		genpar.ibksiz = inp.GetValue("Number of Elements Per Block");

		((string) inp.GetValue("Turn Off Source Terms for Scalars") == "True") ?
				sclrs.nosource = 1 : sclrs.nosource = 0;
		sclrs.tdecay = inp.GetValue("Decay Multiplier for Scalars");

		// TURBULENCE MODELING PARAMETER
		int tpturb = turbvari.iles - turbvari.irans;
		int ifrule;
		if (tpturb != 0) {

			turbvari.nohomog = inp.GetValue("Number of Homogenous Directions");

			if ((string) inp.GetValue("Turbulence Wall Model Type")
					== "Slip Velocity")
				turbvar.itwmod = 1;
			else if ((string) inp.GetValue("Turbulence Wall Model Type")
					== "Effective Viscosity")
				turbvar.itwmod = 2;
			else
				turbvar.itwmod = 0;
			if (turbvari.irans < 0)
				turbvar.itwmod = turbvar.itwmod * (-1);
			ifrule = inp.GetValue("Velocity Averaging Steps");
			turbvar.wtavei = (ifrule > 0) ? 1.0 / ifrule : -1.0 / ifrule;

			if (turbvari.iles == 1) {

				if ((string) inp.GetValue("Dynamic Model Type") == "Bardina")
					turbvari.iles += 10;
				else if ((string) inp.GetValue("Dynamic Model Type")
						== "Projection")
					turbvari.iles += 20;

				ifrule = inp.GetValue("Filter Integration Rule");
				turbvari.iles += ifrule - 1;
				ifrule = inp.GetValue("Dynamic Model Averaging Steps");
				turbvar.dtavei = (ifrule > 0) ? 1.0 / ifrule : -1.0 / ifrule;
				turbvar.fwr1 = inp.GetValue("Filter Width Ratio");
				turbvar.flump = inp.GetValue("Lumping Factor for Filter");

				if ((string) inp.GetValue("Model Statistics") == "True") {
					turbvari.modlstats = 1;
				} else {
					turbvari.modlstats = 0;
				}

				if ((string) inp.GetValue("Double Filter") == "True") {
					turbvari.i2filt = 1;
				} else {
					turbvari.i2filt = 0;
				}

				if ((string) inp.GetValue("Model/SUPG Dissipation") == "True") {
					turbvari.idis = 1;
				} else {
					turbvari.idis = 0;
				}

				if ((string) inp.GetValue("Dynamic Model Type") == "Standard") {

					if ((string) inp.GetValue("Dynamic Sub-Model Type")
							== "None")
						turbvari.isubmod = 0;
					else if ((string) inp.GetValue("Dynamic Sub-Model Type")
							== "DFWR")
						turbvari.isubmod = 1;
					else if ((string) inp.GetValue("Dynamic Sub-Model Type")
							== "SUPG")
						turbvari.isubmod = 2;
				} else if ((string) inp.GetValue("Dynamic Model Type")
						== "Projection") {

					if ((string) inp.GetValue("Projection Filter Type")
							== "Linear")
						turbvari.ifproj = 0;
					else if ((string) inp.GetValue("Projection Filter Type")
							== "Quadratic")
						turbvari.ifproj = 1;

					if ((string) inp.GetValue("Dynamic Sub-Model Type")
							== "None")
						turbvari.isubmod = 0;
					else if ((string) inp.GetValue("Dynamic Sub-Model Type")
							== "ConsistentProj")
						turbvari.isubmod = 1;
				}

			}
		}

		// SPEBC MODELING PARAMETERS

		if ((spebcvr.irscale = inp.GetValue("SPEBC Model Active")) >= 0) {

			ifrule = inp.GetValue("Velocity Averaging Steps");
			turbvar.wtavei =
					(ifrule > 0) ? 1.0 / ifrule : 1.0 / inpdat.nstep[0];
			spebcvr.intpres = inp.GetValue("Interpolate Pressure");
			spebcvr.plandist = inp.GetValue("Distance between Planes");
			spebcvr.thetag = inp.GetValue("Theta Angle of Arc");
			spebcvr.ds = inp.GetValue("Distance for Velocity Averaging");
			spebcvr.tolerence = inp.GetValue("SPEBC Cylindrical Tolerance");
			spebcvr.radcyl = inp.GetValue("Radius of recycle plane");
			spebcvr.rbltin = inp.GetValue("Inlet Boundary Layer Thickness");
			spebcvr.rvscal = inp.GetValue("Vertical Velocity Scale Factor");
		}

		// CARDIOVASCULAR MODELING PARAMETERS
		if ((string) inp.GetValue("Time Varying Boundary Conditions From File")
				== "True")
			nomodule.itvn = 1;
		else
			nomodule.itvn = 0;
		if (nomodule.itvn == 1)
			nomodule.bcttimescale = inp.GetValue("BCT Time Scale Factor");

		nomodule.ipvsq = 0;
		if (nomodule.icardio = inp.GetValue("Number of Coupled Surfaces")) {
			if (nomodule.icardio > MAXSURF) {
				cout << "Number of Coupled Surfaces > MAXSURF \n";
				exit(1);
			}
			if ((string) inp.GetValue("Pressure Coupling") == "None")
				nomodule.ipvsq = 0;
			if ((string) inp.GetValue("Pressure Coupling") == "Explicit")
				nomodule.ipvsq = 1;
			if ((string) inp.GetValue("Pressure Coupling") == "Implicit")
				nomodule.ipvsq = 2;
			if ((string) inp.GetValue("Pressure Coupling") == "P-Implicit")
				nomodule.ipvsq = 3;
			if ((string) inp.GetValue("Inflow Coupling") == "True")
				nomodule.incp = 1;
			else
				nomodule.incp = 0;
			if (nomodule.incp == 1) {
				nomodule.numINCPSrfs = inp.GetValue(
						"Number of Coupled Inflow Surfaces");
				ivec = inp.GetValue("List of Coupled Inflow Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistINCP[i] = 0;
				for (i = 0; i < nomodule.numINCPSrfs; i++) {
					nomodule.nsrflistINCP[i + 1] = ivec[i];
				}
				if ((string) inp.GetValue("Inflow Parameters From File")
						== "True")
					nomodule.incpfile = 1;
				else
					nomodule.incpfile = 0;
			}
			if (nomodule.numResistSrfs = inp.GetValue(
					"Number of Resistance Surfaces")) {
				ivec = inp.GetValue("List of Resistance Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistResist[i] = 0;
				for (i = 0; i < nomodule.numResistSrfs; i++) {
					nomodule.nsrflistResist[i + 1] = ivec[i];
				}
				vec = inp.GetValue("Resistance Values");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.ValueListResist[i] = 0;
				for (i = 0; i < nomodule.numResistSrfs; i++)
					nomodule.ValueListResist[i + 1] = vec[i];
				vec.erase(vec.begin(), vec.end());
			}
			if (nomodule.numImpSrfs = inp.GetValue(
					"Number of Impedance Surfaces")) {
				ivec = inp.GetValue("List of Impedance Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistImp[i] = 0;
				for (i = 0; i < nomodule.numImpSrfs; i++) {
					nomodule.nsrflistImp[i + 1] = ivec[i];
				}
				if ((string) inp.GetValue("Impedance From File") == "True")
					nomodule.impfile = 1;
				else
					nomodule.impfile = 0;
			}
			if (nomodule.numRCRSrfs = inp.GetValue("Number of RCR Surfaces")) {
				ivec = inp.GetValue("List of RCR Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistRCR[i] = 0;
				for (i = 0; i < nomodule.numRCRSrfs; i++) {
					nomodule.nsrflistRCR[i + 1] = ivec[i];
				}
				if ((string) inp.GetValue("RCR Values From File") == "True")
					nomodule.ircrfile = 1;
				else
					nomodule.ircrfile = 0;
			}
			if (nomodule.numCORSrfs = inp.GetValue(
					"Number of Coronary Surfaces")) {
				ivec = inp.GetValue("List of Coronary Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistCOR[i] = 0;
				for (i = 0; i < nomodule.numCORSrfs; i++) {
					nomodule.nsrflistCOR[i + 1] = ivec[i];
				}
				if ((string) inp.GetValue("Coronary Values From File")
						== "True")
					nomodule.icorfile = 1;
				else
					nomodule.icorfile = 0;
			}
			if (nomodule.numVisFluxSrfs = inp.GetValue(
					"Number of Surfaces which zero out in-plane tractions")) {
				ivec = inp.GetValue(
						"List of Surfaces which zero out in-plane tractions");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistVisFlux[i] = 0;
				for (i = 0; i < nomodule.numVisFluxSrfs; i++) {
					nomodule.nsrflistVisFlux[i + 1] = ivec[i];
				}
			}
			if (nomodule.numCalcSrfs = inp.GetValue(
					"Number of Surfaces which Output Pressure and Flow")) {
				ivec = inp.GetValue("List of Output Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistCalc[i] = 0;
				for (i = 0; i < nomodule.numCalcSrfs; i++) {
					nomodule.nsrflistCalc[i + 1] = ivec[i];
				}
			}
			if (nomodule.numDirCalcSrfs =
					inp.GetValue(
							"Number of Dirichlet Surfaces Which Output Pressure and Flow")) {
				ivec = inp.GetValue("List of Dirichlet Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistDirCalc[i] = 0;
				for (i = 0; i < nomodule.numDirCalcSrfs; i++) {
					nomodule.nsrflistDirCalc[i + 1] = ivec[i];
				}
			}
			if ((string) inp.GetValue("Lagrange Multipliers") == "True")
				nomodule.Lagrange = 1;
			else
				nomodule.Lagrange = 0;
			if (nomodule.Lagrange == 1) {
				nomodule.numLagrangeSrfs = inp.GetValue(
						"Number of Constrained Surfaces");
				ivec = inp.GetValue("List of Constrained Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistLagrange[i] = 0;
				for (i = 0; i < nomodule.numLagrangeSrfs; i++) {
					nomodule.nsrflistLagrange[i + 1] = ivec[i];
				}
				if ((string) inp.GetValue(
						"Constrained Surface Information From File") == "True")
					nomodule.iLagfile = 1;
				else
					nomodule.iLagfile = 0;
			}
		}
		nomodule.rescontrol = 0;
		if ((string) inp.GetValue("Residual Control") == "True")
			nomodule.rescontrol = 1;
		if (nomodule.rescontrol == 1) {
			nomodule.ResCriteria = inp.GetValue("Residual Criteria");
			nomodule.MinNumIter = inp.GetValue("Minimum Required Iterations");
		}
		nomodule.ideformwall = 0;
		if ((string) inp.GetValue("Deformable Wall") == "True") {
			nomodule.ideformwall = 1;
			nomodule.rhovw = inp.GetValue("Density of Vessel Wall");
			nomodule.rnuvw = inp.GetValue("Poisson Ratio of Vessel Wall");
			nomodule.rshearconstantvw = inp.GetValue(
					"Shear Constant of Vessel Wall");
			nomodule.nProps = inp.GetValue(
					"Number of Wall Properties per Node");

			if ((string) inp.GetValue("Wall Mass Matrix for LHS") == "True")
				nomodule.iwallmassfactor = 1;
			else
				nomodule.iwallmassfactor = 0;
			if ((string) inp.GetValue("Wall Stiffness Matrix for LHS")
					== "True")
				nomodule.iwallstiffactor = 1;
			else
				nomodule.iwallstiffactor = 0;

			if ((string) inp.GetValue("Use SWB File") == "True")
				nomodule.iUseSWB = 1;
			else {
				nomodule.iUseSWB = 0;
				nomodule.evw = inp.GetValue("Young Mod of Vessel Wall");
				nomodule.thicknessvw = inp.GetValue("Thickness of Vessel Wall");
			}

			if ((string) inp.GetValue("Use TWB File") == "True")
				nomodule.iUseTWB = 1;
			else
				nomodule.iUseTWB = 0;

			if ((string) inp.GetValue("Use EWB File") == "True")
				nomodule.iUseEWB = 1;
			else
				nomodule.iUseEWB = 0;

			if ((string) inp.GetValue("Wall Damping Term") == "True") {
				nomodule.iwalldamp = 1;
				if (nomodule.iUseTWB == 0)
					nomodule.tissSuppDampCoeff = inp.GetValue(
							"Damping Coefficient for Tissue Support");
			} else
				nomodule.iwalldamp = 0;

			if ((string) inp.GetValue("Wall External Support Term") == "True") {
				nomodule.iwallsupp = 1;
				if (nomodule.iUseTWB == 0)
					nomodule.tissSuppStiffCoeff = inp.GetValue(
							"Stiffness Coefficient for Tissue Support");
			} else
				nomodule.iwallsupp = 0;

			if ((string) inp.GetValue("Wall State Filter Term") == "True") {
				nomodule.imeasdist = 1;
				if (nomodule.iUseEWB == 0)
					nomodule.stateFilterCoeff = inp.GetValue(
							"Wall State Filter Coefficient");
			} else
				nomodule.imeasdist = 0;
		}

		// Scaling Parameters Keywords

		outpar.ro = inp.GetValue("Density");
		outpar.vel = inp.GetValue("Velocity");
		outpar.press = inp.GetValue("Pressure");
		outpar.temper = inp.GetValue("Temperature");
		outpar.entrop = inp.GetValue("Entropy");

		// Step Sequencing

		ivec = inp.GetValue("Step Construction");
		sequence.seqsize = ivec.size();
		if (sequence.seqsize > 100 || sequence.seqsize < 2)
			cerr << "Sequence size must be between 2 and 100 " << endl;

		for (i = 0; i < sequence.seqsize; i++)
			sequence.stepseq[i] = ivec[i];

	} catch (exception &e) {
		cout << endl << "Input exception: " << e.what() << endl << endl;
		ierr = 001;
		phSolverMessagePrinter::print_error_code(ierr);
		return ierr;
	}

	return ierr;

}

phSolver::~phSolver() {
}

void phSolverMessagePrinter::print_error_code(int ierr) {
	/*
	 Return Error codes:
	 0xx         Input error
	 1xx         Solution Control
	 105         Turbulence Model not supported

	 2xx         Material Properties

	 3xx         Output Control

	 4xx         Discretization Control

	 5xx         Scaling Parameters

	 6xx         Linear Solver Control
	 601         linear solver type not supported
	 */
	cout << endl << endl << "Input error detected: " << endl << endl;
	if (ierr == 001) {
		cout << endl << "Input Directive not understood" << endl << endl;
	}
	if (ierr == 105) {
		cout << endl << "Turbulence Model Not Supported" << endl << endl;
	}
	if (ierr == 601) {
		cout << endl << "Linear Solver Type Not Supported" << endl << endl;
	}

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
