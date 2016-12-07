#ifndef ABSTRACTBOUNDARYCONDITION_HXX_
#define ABSTRACTBOUNDARYCONDITION_HXX_

#include "gtest/gtest_prod.h"
#include <string>
#include <stdexcept>
#include "FortranBoundaryDataPointerManager.hxx"
#include "debuggingToolsForCpp.hxx"

// Forward declarations:
class BoundaryConditionManager;

class AbstractBoundaryCondition
 {
    friend class testMultidom;
 	FRIEND_TEST(testMultidom, checkBoundaryConditionsMadeProperly);
 	FRIEND_TEST(testMultidom, checkRCRLinearInterpolators);
 	// FRIEND_TEST(testMultidom, checkDpDqAndHopFortranPasser)
 	FRIEND_TEST(testMultidom, checkImplicitConditionComputation_solve);
 	FRIEND_TEST(testMultidom, checkImplicitConditionComputation_update);
 	FRIEND_TEST(testMultidom, checkFlowAndPressureSetters);
 public:
    AbstractBoundaryCondition(const int surfaceIndex_in, const double hstep_in, const double delt_in, const double alfi_in, const double currentTimestepIndex, const int maxsurf, const int nstep)
    : surfaceIndex(surfaceIndex_in),
      hstep(hstep_in),
      delt(delt_in),
      m_generalizedAlphaMethodAlpha(alfi_in),
      m_currentTimestepIndex(currentTimestepIndex),
      m_maxsurf(maxsurf),
      m_nstep(nstep)
    {
        hasListOfMeshNodesAtThisBoundary = false; // flag to be used to guard against using the list when it hasn't been provided by Fortran.
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

    virtual ~AbstractBoundaryCondition()
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
    double getSurfaceIndex() const {return surfaceIndex;}
    void incrementTimestepIndex();
    virtual void setDirichletConditionsIfNecessary(int* const binaryMask);
    virtual void resetStateUsingKalmanFilteredEstimate(const double flow, const double pressure, const int timestepNumber) = 0;
    void debugPrintFlowPointerTarget();
    void setFlowN(const double flowN);
    void setFlowN1(const double flowN);
    void setListOfMeshNodesAtThisBoundary(const int* const & ndsurf_nodeToBoundaryAssociationArray, const int& lengthOfNodeToBoundaryAssociationArray);
    void updatePressureAndFlowHistory();
    double getPressureHistoryValueByTimestepIndex(const int timestepIndex) const;
    double getFlowHistoryValueByTimestepIndex(const int timestepIndex) const;
    void computeImplicitCoeff_solve(const int timestepNumber);
    void computeImplicitCoeff_update(const int timestepNumber);
    virtual void finaliseAtEndOfTimestep() = 0;
 protected:
 	double dp_dq;
 	double Hop;
 	double dp_dq_n1;
 	double Hop_n1;
 	double pressure_n;
 	double flow_n;
 	const int surfaceIndex;
    const int hstep;
    const double delt;
    const double m_generalizedAlphaMethodAlpha;
    int m_currentTimestepIndex;
    const int m_maxsurf;
    const int m_nstep;
 	int isactive;
 	double* flowhist;
 	double* pressurehist;
 	std::string flowfile;
    std::string pressurefile;
	std::vector<double*> flow_n_ptrs;
    double flow_n1;
    std::vector<double*> pressure_n_ptrs;
    // double implicitcoeff;
    // double implicitcoeff_n1; 

    bool m_thisIsARestartedSimulation;
    std::vector<int> listOfMeshNodesAtThisBoundary;
    bool hasListOfMeshNodesAtThisBoundary;

    // double LPNInflowPressure;

 private:
 	static int numberOfConstructedBoundaryConditions;
    virtual std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt) = 0;
 };

#endif