#ifndef ABSTRACTBOUNDARYCONDITION_HXX_
#define ABSTRACTBOUNDARYCONDITION_HXX_

#include "gtest/gtest_prod.h"
#include <string>
#include <stdexcept>
#include "fortranPointerManager.hxx"
#include "debuggingToolsForCpp.hxx"

// Forward declarations:
class boundaryConditionManager;

class abstractBoundaryCondition
 {
 	friend class boundaryConditionManager;
 	friend class testMultidom;
 	friend class basicFileWriter;
 	FRIEND_TEST(testMultidom, checkBoundaryConditionsMadeProperly);
 	FRIEND_TEST(testMultidom, checkRCRLinearInterpolators);
 	// FRIEND_TEST(testMultidom, checkDpDqAndHopFortranPasser)
 	FRIEND_TEST(testMultidom, checkImplicitConditionComputation_solve);
 	FRIEND_TEST(testMultidom, checkImplicitConditionComputation_update);
 	FRIEND_TEST(testMultidom, checkFlowAndPressureSetters);
 public:
    abstractBoundaryCondition(const int surfaceIndex_in, const double hstep_in, const double delt_in, const double alfi_in, const double currentTimestepIndex, const int maxsurf, const int nstep)
    : surfaceIndex(surfaceIndex_in),
      hstep(hstep_in),
      delt(delt_in),
      alfi_local(alfi_in),
      m_currentTimestepIndex(currentTimestepIndex),
      m_maxsurf(maxsurf),
      m_nstep(nstep)
    {
        hasListOfMeshNodesAtThisBoundary = false; // flag to be used to guard against using the list when it hasn't been provided by Fortran.
        std::cout <<"is surfarea set yet in c++?" << surfarea << std::endl;
        // allocate arrays with +1 to size, in case hstep=0 (that would be undefined behaviour under new double)
        flowhist = new double [hstep+1];
        pressurehist = new double [hstep+1];

        dp_dq = 0.0;
        Hop = 0.0;
        dp_dq_n1 = 0.0;
        Hop_n1 = 0.0;
        numberOfConstructedBoundaryConditions++;
        index = numberOfConstructedBoundaryConditions;

        if (m_currentTimestepIndex > 0)
        {
            m_thisIsARestartedSimulation = true;
        }
        else
        {
            m_thisIsARestartedSimulation = false;
        }

    }

    virtual void getPressureAndFlowPointersFromFortran();
    virtual bool flowPermittedAcross3DInterface();

    virtual ~abstractBoundaryCondition()
    {
        delete[] flowhist;
        delete[] pressurehist;
        numberOfConstructedBoundaryConditions--;
    }
    virtual void initialiseModel() = 0;
    // virtual void initialiseAtStartOfTimestep() = 0;
    double getdp_dq();
    double getHop();
    int index;

    // void setLPNInflowPressure(double inflowPressure);
    void updpressure_n1_withflow();
    // virtual void finalizeLPNAtEndOfTimestep() = 0;
    double getSurfaceIndex() const {return surfaceIndex;}
    void incrementTimestepIndex();
    virtual void setDirichletConditionsIfNecessary(int* const binaryMask);
    // void storeFlowAndPressureAtStartOfTimestep();
 protected:
 	double dp_dq;
 	double Hop;
 	double dp_dq_n1;
 	double Hop_n1;
 	double pressure_n;
 	double flow_n;
 	int surfaceIndex;
 	int isactive;
 	double* flowhist;
 	double* pressurehist;
 	std::string flowfile;
    std::string pressurefile;
	double surfarea;
	std::vector<double*> flow_n_ptrs;
    double flow_n1;
    std::vector<double*> pressure_n_ptrs;
    // double implicitcoeff;
    // double implicitcoeff_n1; 
    const int hstep;
    int m_currentTimestepIndex;
    const double delt;
    const double alfi_local;
    const int m_maxsurf;
    const int m_nstep;
    bool m_thisIsARestartedSimulation;
    std::vector<int> listOfMeshNodesAtThisBoundary;
    bool hasListOfMeshNodesAtThisBoundary;

    // double LPNInflowPressure;

    void setListOfMeshNodesAtThisBoundary(const int* const & ndsurf_nodeToBoundaryAssociationArray, const int& lengthOfNodeToBoundaryAssociationArray);
    void computeImplicitCoeff_solve(const int timestepNumber);
 	void computeImplicitCoeff_update(const int timestepNumber);
    virtual std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt) = 0;

	void updatePressureAndFlowHistory();
    // virtual void updateLPN() = 0;
	virtual double linInterpolateTimeData(const double &currentTime, const int timeDataLength)
	{
    	throw std::runtime_error("Disallowed call to non-overridden (e.g. non-RCR).");
    	return 0.0;
    };
 private:
 	static int numberOfConstructedBoundaryConditions;
 };

#endif