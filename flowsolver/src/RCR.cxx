#include "RCR.hxx"
#include "FortranBoundaryDataPointerManager.hxx"
#include "BoundaryConditionManager.hxx"

// Statics
int RCR::s_numberOfInitialisedRCRs = 0;

/*  >>>>>>>>>>>>>>>>>>>>     RUUUULE BRITANNIAAAA! <<<<<<<<<<<<<<<<<<<<<<<<<<
ZZZ         888888888888888888888   ZZZZZZZZ  888888888888888888888   ZZZZZZ    
ZZZZZZ         888888888888888888   ZZZZZZZZ  888888888888888888   ZZZZZZ       
   ZZZZZZ         888888888888888   ZZZZZZZZ  888888888888888   ZZZZZZ         8
888   ZZZZZZ         888888888888   ZZZZZZZZ  888888888888   ZZZZZZ         8888
8888887  =ZZZZZ:        888888888   ZZZZZZZZ  888888888   ZZZZZZ        :8888888
8888888888   ZZZZZZ         88888   ZZZZZZZZ  88888   ZZZZZZ         88888888888
8888888888888   ZZZZZZ         88   ZZZZZZZZ  88   ZZZZZZ         88888888888888
8888888888888888   ZZZZZZ           ZZZZZZZZ    ZZZZZZ         88888888888888888
                                    ZZZZZZZZ                                    
                                    ZZZZZZZZ                                    
ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
                                    ZZZZZZZZ                                    
8888888888888888         ZZZZZZ     ZZZZZZZZ          ZZZZZZ   88888888888888888
8888888888888         ZZZZZZ   88   ZZZZZZZZ  88         ZZZZZZ   88888888888888
8888888888         ZZZZZZ   88888   ZZZZZZZZ  88888         ZZZZZZ   88888888888
8888888         ZZZZZ$  :88888888   ZZZZZZZZ  888888887        IZZZZZ.  88888888
888         ZZZZZZ   888888888888   ZZZZZZZZ  888888888888         ZZZZZZ   8888
         ZZZZZZ   888888888888888   ZZZZZZZZ  888888888888888         ZZZZZZ   8
      ZZZZZZ   888888888888888888   ZZZZZZZZ  888888888888888888         ZZZZZZ 
   ZZZZZZ   888888888888888888888   ZZZZZZZZ  888888888888888888888         ZZZZ
ZZZZZ$  :888888888888888888888888   ZZZZZZZZ  8888888888888888888888887        I
*/
void RCR::initialiseModel()
{
    // History arrays      
    for(int ii=0; ii<hstep; ii++)
    {
        flowhist[ii] = 0.0;
        pressurehist[ii] = 0.0;
    }

    flowfile = "QHistRCR.dat";
    pressurefile = "PHistRCR.dat";
}

// Here we step the actual discretised ODE for the RCR:
std::pair<double,double> RCR::computeImplicitCoefficients(const int timestepNumber, const double timeAtStepNplus1, const double alfi_delt)
{
  double timeAtStepN = delt*((double)timestepNumber);

  double pdistn = m_pDistLinearInterpolator.interpolateInTimeWithConstantExtrapolation(timeAtStepN);
  double pdistn_1 = m_pDistLinearInterpolator.interpolateInTimeWithConstantExtrapolation(timeAtStepNplus1);

  // // parameters overwritten
  // // dirty hack for filtering
  //  distalResistance = a%parameters_RCR(3,i)
  //  proximalResistance = a%parameters_RCR(1,i)
  //  capacitance = a%parameters_RCR(2,i)

  //  pdistn = a%parameters_Pd
  //  pdistn_1 = a%parameters_Pd

  // std::cout << "C++ RCR compliance: " << capacitance << " rp: "  << proximalResistance << " rd: " << distalResistance << " timestep: " << timestepNumber << " flow: " << (*flow_n_ptrs.at(0))<< " pressure: " << pressure_n << std::endl;
  // std::cout << "pdistn: " << pdistn << "pdistn_1: " << pdistn_1 << std::endl;

  std::pair<double,double> returnCoeffs;

  double denominator = 1.0 + ((capacitance*distalResistance)/alfi_delt);

  double firstCoefficientNumerator = distalResistance + proximalResistance * (1.0 + ((capacitance*distalResistance)/alfi_delt));
  returnCoeffs.first = firstCoefficientNumerator / denominator;

  double temp2 = pressure_n + pdistn_1 - pdistn - proximalResistance * (*flow_n_ptrs.at(0));
  temp2 = ((capacitance*distalResistance)/alfi_delt)*temp2+ pdistn_1;

  returnCoeffs.second = temp2 / denominator;

  return returnCoeffs;
}

// This function exists because I can't untangle the itrdrv_init in itrdrv.f90 in
// such a way that this gets set up naturally. The Fortran side of things needs
// some major refactoring, but for now, we have to make do with this hack to
// initialise the pressure.
void RCR::setPressureFromFortran()
{
  // This is only called if this is a new simulation (the bool sees to that).
  // If it's a restarted sim, the pressure should be loaded properly anyway.
  if (m_needsPressureToBeInitialisedFromFortran)  // possibly need to disable this guard for kalman filter
  {
    pressure_n = *(pressure_n_ptrs.at(0));
  }
}

void RCR::resetStateUsingKalmanFilteredEstimate(const double flow, const double pressure, const int timestepNumber)
{
  std::cout << "rcr setting timestepNumber: " << timestepNumber << std::endl;
  if (timestepNumber > 0) {
    // These values are used to specify the /history/ state of the model (i.e. P_{n-1} and Q_{n-1})
    *flow_n_ptrs.at(0) = flow;
    pressure_n = pressure;
  }
}

// This overrides the base class version of this function call; on a restart, for the RCR, we need to make sure
// we don't overwrite the pressure data from the pressure history file with an irrelevant value from Fortran.
void RCR::getPressureAndFlowPointersFromFortran()
{
    // here we set the initial values of the flow and pressure using the pointers to the multidomaincontainer.
    // NB: Need to add a method in fortran to set a value for non-zero restarting!
    assert(flow_n_ptrs.size()==0);
    double* flowPointer = FortranBoundaryDataPointerManager::Get()->getBoundaryFlows(surfaceIndex);
    flow_n_ptrs.push_back(flowPointer);
    assert(pressure_n_ptrs.size()==0);
    pressure_n_ptrs.push_back(FortranBoundaryDataPointerManager::Get()->getBoundaryPressures(surfaceIndex));

    flow_n = *(flow_n_ptrs.at(0));
    flow_n1 = 0.0;

    if (m_thisIsARestartedSimulation)
    {
      // Initialise the pressure using the value from the PHistRCR.dat.
      pressure_n = (BoundaryConditionManager::Instance()->PHistReader)->getReadFileData(m_indexOfThisRCR+1, m_currentTimestepIndex);
    }
    else
    {
      m_needsPressureToBeInitialisedFromFortran = true;
    }
}