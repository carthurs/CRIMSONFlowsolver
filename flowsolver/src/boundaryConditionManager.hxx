#ifndef BOOUNDARYCONDITIONMANAGER_HXX_
#define BOOUNDARYCONDITIONMANAGER_HXX_

#include "fileReaders.hxx"
#include "boundaryConditionFactory.hxx"
#include "abstractBoundaryCondition.hxx"
#include "ControlSystemsManager.hxx"
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
	
    void setSurfaceList(const std::vector<std::pair<int,std::string>> surfaceList);

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

    void computeAllImplicitCoeff_solve(const int timestepNumber);
    void computeAllImplicitCoeff_update(const int timestepNumber);

    void updateAllRCRS_Pressure_n1_withflow();
    void updateAllRCRS_setflow_n(const double* const flows);
    void updateAllRCRS_setflow_n1(const double* const flows);

    void recordPressuresAndFlowsInHistoryArrays();

    void writePHistAndQHistRCR();

    // void setSurfacePressure_controlledCoronary(double* coronarySurfacePressures);
    void getImplicitCoeff_controlledCoronary(double* const implicitCoeffs_toBeFilled);
    void updateAllControlledCoronaryLPNs();
    void finalizeLPNAtEndOfTimestep_controlledCoronary();
    void finalizeLPNAtEndOfTimestep_netlists();

    // void updateAllControlledCoronaryLPNs_Pressure_n1_withflow();

    // void setSurfacePressure_netlistLPNs(double* netlistSurfacePressures);
    void initialiseLPNAtStartOfTimestep_netlist();
    void updateAllNetlistLPNs(const int timestepNumber);
    void getImplicitCoeff_netlistLPNs(double* const implicitCoeffs_toBeFilled);
    std::map<int,std::pair<double,double>> getImplicitCoeff_netlistLPNs_toPassTo3DDomainReplacement();
    void writeAllNetlistComponentFlowsAndNodalPressures();

    int getNumberOfRCRSurfaces(){return m_NumberOfRCRSurfaces;}
    int getNumberOfNetlistSurfaces(){return m_NumberOfNetlistSurfaces;}
    int getNumberOfControlledCoronarySurfaces(){return m_NumberOfControlledCoronarySurfaces;}

    void setNumberOfRCRSurfaces(const int numGRCRSrfs);
    void setNumberOfControlledCoronarySurfaces(const int numControlledCoronarySrfs);
    void setNumberOfNetlistSurfaces(const int numNetlistLPNSrfs);
    void setDelt(const double delt);
    void setHstep(const int hstep);
    void setAlfi(const double alfi);
    void setLstep(const int lstep);
    void setNtout(const int ntout);
    void setMaxsurf(const int maxsurf);
    void setNstep(const int nstep);

    void createControlSystems();
    void updateAllControlSystems();

    std::vector<std::pair<boundary_data_t,double>> getBoundaryPressuresOrFlows_zeroDDomainReplacement(const int timestepNumber);

    void setZeroDDomainReplacementPressuresAndFlows(double* zeroDDomainPressures, double* zeroDDomainFlows);

    ~boundaryConditionManager()
    {
    }
 
 private:
    boundaryConditionManager()
    {
        checkIfThisIsARestartedSimulation();
        m_nextTimestepWrite_netlistBoundaries_start = 0;

        m_deltHasBeenSet = false;
        m_hstepHasBeenSet = false;
        m_alfiHasBeenSet = false;
        m_lstepHasBeenSet = false;
        m_maxsurfHasBeenSet = false;
        m_nstepHasBeenSet = false;
    }
    std::vector<boost::shared_ptr<abstractBoundaryCondition>> boundaryConditions;
    // std::map<int,std::pair<double,double>> implicitCoefficientMap;

    std::unique_ptr<ControlSystemsManager> mp_controlSystemsManager;

    static bool m_thisIsARestartedSimulation;
    int m_NumberOfRCRSurfaces;
    int m_NumberOfNetlistSurfaces;
    int m_NumberOfControlledCoronarySurfaces;
    double m_delt;
    int m_hstep;
    double m_alfi;
    int m_lstep;
    int m_ntout;
    int m_maxsurf;
    int m_nstep;

    bool m_deltHasBeenSet;
    bool m_hstepHasBeenSet;
    bool m_alfiHasBeenSet;
    bool m_lstepHasBeenSet;
    bool m_ntoutHasBeenSet;
    bool m_maxsurfHasBeenSet;
    bool m_nstepHasBeenSet;

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