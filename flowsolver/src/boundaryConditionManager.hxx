#ifndef BOOUNDARYCONDITIONMANAGER_HXX_
#define BOOUNDARYCONDITIONMANAGER_HXX_

#include "fileReaders.hxx"
#include "boundaryConditionFactory.hxx"
#include "abstractBoundaryCondition.hxx"
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
    		instance->ifRestartingLoadNecessaryData();
		}
		return instance;
	}

	static void Term()
	{
		if (thisIsARestartedSimulation)
		{
			delete boundaryConditionManager::PHistReader;
		}
		delete instance;
		instance = 0;
	}
	
    void setSurfaceList(const std::vector<std::pair<int,std::string>> surfaceList);

    void ifRestartingLoadNecessaryData();
    
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

    // void updateAllControlledCoronaryLPNs_Pressure_n1_withflow();

    // void setSurfacePressure_netlistLPNs(double* netlistSurfacePressures);
    void initialiseLPNAtStartOfTimestep_netlist();
    void updateAllNetlistLPNs();
    void getImplicitCoeff_netlistLPNs(double* const implicitCoeffs_toBeFilled);
    void writeAllNetlistComponentFlowsAndNodalPressures();

    int getNumberOfRCRSurfaces(){return numberOfRCRSurfaces;}
    int getNumberOfNetlistSurfaces(){return numberOfNetlistSurfaces;}
    int getNumberOfControlledCoronarySurfaces(){return numberOfControlledCoronarySurfaces;}

    ~boundaryConditionManager()
    {
    }
 
 private:
    boundaryConditionManager()
    {
    	numberOfRCRSurfaces = grcrbccom.numGRCRSrfs;
        numberOfNetlistSurfaces = nomodule.numNetlistLPNSrfs;
        numberOfControlledCoronarySurfaces = nomodule.numControlledCoronarySrfs;
    	if (timdat.lstep > 0)
    	{
    		thisIsARestartedSimulation = 1;
            SimpleFileReader numstartReader("numstart.dat");

            bool success = false;
            std::string numstartString = numstartReader.getNextDataSplitBySpacesOrEndOfLine(success);
            assert(success);

            m_nextTimestepWrite_start = boost::lexical_cast<int>(numstartString)+1; // +1 because numstart should contain the step just written before the program last terminated. So we need to start writing on the next (+1 th) time-step.
    	}
    	else
    	{
    		thisIsARestartedSimulation = 0;
            m_nextTimestepWrite_start = 0;
    	}
    }
    std::vector<boost::shared_ptr<abstractBoundaryCondition>> boundaryConditions;
    // std::map<int,std::pair<double,double>> implicitCoefficientMap;
    boundaryConditionFactory factory;

    static int thisIsARestartedSimulation;
    int numberOfRCRSurfaces;
    int numberOfNetlistSurfaces;
    int numberOfControlledCoronarySurfaces;

    int m_nextTimestepWrite_start;
    int m_nextTimestepWrite_end;

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