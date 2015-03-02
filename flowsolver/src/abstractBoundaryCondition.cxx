#include "abstractBoundaryCondition.hxx"
#include <stdexcept>


// Statics
int abstractBoundaryCondition::numberOfConstructedBoundaryConditions = 0;

// initialise the multidomain/LPN objects, this will need IFDEF for 3D and 1D codes
double abstractBoundaryCondition::getHop()
{
  return Hop;
}

double abstractBoundaryCondition::getdp_dq()
{
  return dp_dq;
}

// void abstractBoundaryCondition::setLPNInflowPressure(double inflowPressure)
// {
// 	// This will not be needed by all subclasses!
//     LPNInflowPressure = inflowPressure;
// }

void abstractBoundaryCondition::computeImplicitCoeff_solve(const int timestepNumber)
{
  std::pair<double,double> temp;

  double timeAtStepNplus1 = delt*((double)timestepNumber+alfi_local);
  double alfi_delt = alfi_local*delt;

  temp = computeImplicitCoefficients(timestepNumber, timeAtStepNplus1, alfi_delt);

  dp_dq = temp.first;
  Hop = temp.second;
}

void abstractBoundaryCondition::computeImplicitCoeff_update(const int timestepNumber)
{
  std::pair<double,double> temp;

  double timeAtStepNplus1 = delt*((double)timestepNumber+1.0);
  double alfi_delt = delt;

  temp = computeImplicitCoefficients(timestepNumber, timeAtStepNplus1, alfi_delt);

  dp_dq_n1 = temp.first;
  Hop_n1 = temp.second;
}

void abstractBoundaryCondition::updatePressureAndFlowHistory()
{
  pressurehist[timdat.lstep] = pressure_n;
  flowhist[timdat.lstep] = *flow_n_ptr;
}

void abstractBoundaryCondition::initialiseAtStartOfTimestep()
{
  throw std::logic_error("EE: Call to initialiseAtStartOfTimestep as a member of abstractBoundaryCondition. This is not allowed.");
}

void abstractBoundaryCondition::updateLPN()
{
	throw std::logic_error("Error: Call to updateLPN as a member of abstractBoundaryCondition. This is not allowed.");
}

void abstractBoundaryCondition::updpressure_n1_withflow()
{
  pressure_n = dp_dq_n1*(*flow_n_ptr) + Hop_n1;
}

void abstractBoundaryCondition::finalizeLPNAtEndOfTimestep()
{
	throw std::logic_error("Error: Call to updateLPN as a member of abstractBoundaryCondition. This is not allowed.");
}

void abstractBoundaryCondition::setListOfMeshNodesAtThisBoundary(const int* const & ndsurf_nodeToBoundaryAssociationArray, const int& lengthOfNodeToBoundaryAssociationArray)
{
  for (int node=0; node<lengthOfNodeToBoundaryAssociationArray; node++)
  {
    if (ndsurf_nodeToBoundaryAssociationArray[node] == surfaceIndex)
    {
      listOfMeshNodesAtThisBoundary.push_back(node);
    }
  }
  hasListOfMeshNodesAtThisBoundary = true;   
}

void abstractBoundaryCondition::getPressureAndFlowPointersFromFortran()
{
    // here we set the initial values of the flow and pressure using the pointers to the multidomaincontainer.
    // NB: Need to add a method in fortran to set a value for non-zero restarting!
    flow_n_ptr = fortranBoundaryDataPointerManager::Get()->getBoundaryFlows(surfaceIndex);
    pressure_n_ptr = fortranBoundaryDataPointerManager::Get()->getBoundaryPressures(surfaceIndex);

    flow_n = *flow_n_ptr;
    flow_n1 = 0.0;
    pressure_n = *pressure_n_ptr;
}