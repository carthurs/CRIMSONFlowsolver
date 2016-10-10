#ifndef BOOUNDARYCONDITIONMANAGER_HXX_
#define BOOUNDARYCONDITIONMANAGER_HXX_

#include "fileReaders.hxx"
#include "boundaryConditionFactory.hxx"
#include "abstractBoundaryCondition.hxx"
#include "ControlSystemsManager.hxx"
#include "ClosedLoopDownstreamSubsection.hxx"
#include <boost/lexical_cast.hpp>

// Forward declarations:
class abstractBoundaryCondition;

 class boundaryConditionManager
 {
 	friend class testMultidom;
 public:

 	static histFileReader* PHistReader;
 
    static boundaryConditionManager* Instance()
	{
		if (!instance) {
			instance = new boundaryConditionManager();
		}
		return instance;
	}

	static void Term()
	{
		if (m_thisIsARestartedSimulation)
		{
            if (boundaryConditionManager::PHistReader != NULL)
            {
			     delete boundaryConditionManager::PHistReader;
                 boundaryConditionManager::PHistReader = NULL;
            }
		}
		delete instance;
		instance = 0;
	}
	
    void setSurfaceList(const std::vector<std::pair<int,boundary_condition_t>> surfaceList);

    void ifRestartingLoadNecessaryData();

    void giveBoundaryConditionsListsOfTheirAssociatedMeshNodes(const int* ndsurf_nodeToBoundaryAssociationArray, const int& lengthOfNodeToBoundaryAssociationArray);
    void getBinaryMaskToAdjustNodalBoundaryConditions(int* const binaryMask, const int binaryMaskLength);
    void getNumberOfBoundaryConditionsWhichCurrentlyDisallowFlow(int& numBCsWhichDisallowFlow);
    void getNumberOfNetlistBoundaryConditionsWhichCurrentlyAllowFlow(int& numBCsWhichDisallowFlow);
    void discoverWhetherFlowPermittedAcrossSurface(const int& queriedSurfaceIndex, int& flowIsPermitted);
    void haveBoundaryConditionTypesChanged(int& boundaryConditionTypesHaveChanged);
    
    void setPressureFromFortran();
    void getImplicitCoeff_rcr(double* const implicitCoeffs_toBeFilled);
    std::vector<boost::shared_ptr<abstractBoundaryCondition>>* getBoundaryConditions();

    // Must be fully defined here in the header so that other translation units can instantiate the template
    template <typename TemplateBoundaryConditionType>
    void computeImplicitCoeff_solve(const int timestepNumber)
    {
      for (auto&& boundaryCondition : m_boundaryConditions)
      {
        if (boost::dynamic_pointer_cast<TemplateBoundaryConditionType> (boundaryCondition))
        {
          boundaryCondition->computeImplicitCoeff_solve(timestepNumber);
        }
      }
    }

    template <typename TemplateBoundaryConditionType>
    void computeImplicitCoeff_update(const int timestepNumber);

    void updateAllRCRS_Pressure_n1_withflow();
    void updateAllRCRS_setflow_n(const double* const flows);
    void updateAllRCRS_setflow_n1(const double* const flows);

    // void storeAllBoundaryConditionFlowsAndPressuresAtStartOfTimestep();

    void recordPressuresAndFlowsInHistoryArrays();

    void writePHistAndQHistRCR();

    // void setSurfacePressure_controlledCoronary(double* coronarySurfacePressures);
    void getImplicitCoeff_controlledCoronary(double* const implicitCoeffs_toBeFilled);
    void updateAllControlledCoronaryLPNs();
    void finalizeLPNAtEndOfTimestep_controlledCoronary();
    void finalizeLPNAtEndOfTimestep_netlists();

    std::vector<double*> getPointersToAllNetlistCapacitorNodalHistoryPressures() const;

    // void updateAllControlledCoronaryLPNs_Pressure_n1_withflow();

    // void setSurfacePressure_netlistLPNs(double* netlistSurfacePressures);
    void initialiseLPNAtStartOfTimestep_netlist();
    void updateAllNetlistLPNs(const int timestepNumber);
    void getImplicitCoeff_netlistLPNs(double* const implicitCoeffs_toBeFilled);
    std::map<int,std::pair<double,double>> getImplicitCoeff_netlistLPNs_toPassTo3DDomainReplacement();
    void writeAllNetlistComponentFlowsAndNodalPressures();
    // void loadAllNetlistComponentFlowsAndNodalPressures();

    int getNumberOfRCRSurfaces(){return m_NumberOfRCRSurfaces;}
    int getNumberOfNetlistSurfaces(){return m_NumberOfNetlistSurfaces;}
    int getNumberOfControlledCoronarySurfaces(){return m_NumberOfControlledCoronarySurfaces;}
    int getNumberOfDownsreamClosedLoopCircuits(){return m_numLoopClosingNetlistCircuits;}
    void getNumberOfBoundaryConditionManagerBoundaryConditions_reference(int& totalNumberOfManagedBoundaryConditions) const;

    void setNumberOfRCRSurfaces(const int numGRCRSrfs);
    void setNumberOfControlledCoronarySurfaces(const int numControlledCoronarySrfs);
    void setNumberOfNetlistSurfaces(const int numNetlistLPNSrfs);
    void setDelt(const double delt);
    void setHstep(const int hstep);
    void setAlfi(const double alfi);
    // void setLstep(const int currentTimestepIndex);
    void setStartingTimestepIndex(const int startingTimestepIndex);
    void incrementTimestepIndex();
    void setNtout(const int ntout);
    void setMaxsurf(const int maxsurf);
    void setNstep(const int nstep);
    void setNumLoopClosingnetlistCircuits(const int numLoopClosingCircuits);
    void setMasterControlScriptPresent(const int masterControlScriptPresent);

    void resetStateUsingKalmanFilteredEstimate(const double flow, const double pressure, const int surfaceIndex, const int timestepNumber);

    void setSimulationModePurelyZeroD(const int simulationIsPurelyZeroD);

    void createControlSystems();
    void updateBoundaryConditionControlSystems();

    void markClosedLoopLinearSystemsForRebuilding();

    std::vector<std::pair<boundary_data_t,double>> getBoundaryPressuresOrFlows_zeroDDomainReplacement(const int timestepNumber);

    void setZeroDDomainReplacementPressuresAndFlows(double* zeroDDomainPressures, double* zeroDDomainFlows);

    void debugPrintFlowPointerTarget_BCM();

    ~boundaryConditionManager()
    {
    }
 
 private:
    boundaryConditionManager()
    {
        m_numberOfBoundaryConditionsManaged = 0;

        m_NumberOfRCRSurfaces = 0;
        m_NumberOfNetlistSurfaces = 0;
        m_NumberOfControlledCoronarySurfaces = 0;

        m_deltHasBeenSet = false;
        m_hstepHasBeenSet = false;
        m_alfiHasBeenSet = false;
        m_currentTimestepIndexHasBeenSet = false;
        m_maxsurfHasBeenSet = false;
        m_nstepHasBeenSet = false;
        m_numLoopClosingNetlistCircuitsHasBeenSet = false;
        m_startingTimestepIndexHasBeenSet = false;

        m_hasSurfaceList = false;
        m_controlSystemsPresent = false;
        m_masterControlScriptPresent = false;

        m_simulationIsPurelyZeroD = false;

        checkIfThisIsARestartedSimulation();
    }
    std::vector<boost::shared_ptr<abstractBoundaryCondition>> m_boundaryConditions;
    std::vector<boost::shared_ptr<ClosedLoopDownstreamSubsection>> m_netlistDownstreamLoopClosingSubsections;
    // std::map<int,std::pair<double,double>> implicitCoefficientMap;

    boost::shared_ptr<ControlSystemsManager> mp_controlSystemsManager;

    bool m_simulationIsPurelyZeroD;

    static bool m_thisIsARestartedSimulation;
    int m_numberOfBoundaryConditionsManaged;
    int m_NumberOfRCRSurfaces;
    int m_NumberOfNetlistSurfaces;
    int m_NumberOfControlledCoronarySurfaces;
    double m_delt;
    int m_hstep;
    double m_alfi;
    int m_currentTimestepIndex; // Formerly "currentTimestepIndex", but that was a silly name.
    int m_ntout;
    int m_maxsurf;
    int m_nstep;
    int m_startingTimestepIndex;
    int m_numLoopClosingNetlistCircuits;

    bool m_deltHasBeenSet;
    bool m_hstepHasBeenSet;
    bool m_alfiHasBeenSet;
    bool m_currentTimestepIndexHasBeenSet;
    bool m_ntoutHasBeenSet;
    bool m_maxsurfHasBeenSet;
    bool m_nstepHasBeenSet;
    bool m_numLoopClosingNetlistCircuitsHasBeenSet;
    bool m_masterControlScriptPresent;
    bool m_startingTimestepIndexHasBeenSet;

    bool m_hasSurfaceList;
    bool m_controlSystemsPresent;

    int m_nextTimestepWrite_netlistBoundaries_start;

    void checkIfThisIsARestartedSimulation();

    // For testing purposes, to clear the static class out before the next test begins
    // Note that you'll have to make the test class a friend in order to use this..
    // I've made it private on purpose!
    void tearDown()
    {
		boundaryConditionManager::Instance()->Term();
    }

    static boundaryConditionManager* instance;
 };

 #endif