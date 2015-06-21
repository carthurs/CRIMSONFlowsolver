#include "abstractBoundaryCondition.hxx"
#include <stdexcept>
#include <cmath>


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

bool abstractBoundaryCondition::flowPermittedAcross3DInterface()
{
  // This is only ever false in some cases, where some subclass of abstractBoundaryCondition overrides this method.
  return true;
}

void abstractBoundaryCondition::computeImplicitCoeff_solve(const int timestepNumber)
{
  if (flowPermittedAcross3DInterface())
  {
    std::pair<double,double> temp;

    double timeAtStepNplus1 = delt*((double)timestepNumber+alfi_local);
    double alfi_delt = alfi_local*delt;

    temp = computeImplicitCoefficients(timestepNumber, timeAtStepNplus1, alfi_delt);

    dp_dq = temp.first;
    Hop = temp.second;
  }
  else
  {
    // else this is a case where the boundary conditions are Dirichlet, so the Hop and dp_dq will not be used.
    // Set them to NaN, defensively.
    dp_dq = NAN;
    Hop = NAN;
  }
}

void abstractBoundaryCondition::computeImplicitCoeff_update(const int timestepNumber)
{
  if (flowPermittedAcross3DInterface())
  {
    std::pair<double,double> temp;

    double timeAtStepNplus1 = delt*((double)timestepNumber+1.0);

    temp = computeImplicitCoefficients(timestepNumber, timeAtStepNplus1, delt);

    dp_dq_n1 = temp.first;
    Hop_n1 = temp.second;
  }
  else
  {
    // else this is a case where the boundary conditions are Dirichlet, so the Hop and dp_dq will not be used.
    // Set them to NaN, defensively.
    dp_dq_n1 = NAN;
    Hop_n1 = NAN;
  }
}

void abstractBoundaryCondition::updatePressureAndFlowHistory()
{
  pressurehist[m_currentTimestepIndex] = pressure_n;
  flowhist[m_currentTimestepIndex] = *(flow_n_ptrs.at(0));
}

void abstractBoundaryCondition::updpressure_n1_withflow()
{
  pressure_n = dp_dq_n1*(*(flow_n_ptrs.at(0))) + Hop_n1;
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
    assert(flow_n_ptrs.size()==0);
    flow_n_ptrs.push_back(fortranBoundaryDataPointerManager::Get()->getBoundaryFlows(surfaceIndex));
    assert(pressure_n_ptrs.size()==0);
    pressure_n_ptrs.push_back(fortranBoundaryDataPointerManager::Get()->getBoundaryPressures(surfaceIndex));

    flow_n = *(flow_n_ptrs.at(0));
    flow_n1 = 0.0;
    pressure_n = *(pressure_n_ptrs.at(0));
}

void abstractBoundaryCondition::incrementTimestepIndex()
{
  m_currentTimestepIndex++;
}

// Processes the binaryMask for setting Dirichlet conditions.
// This boundary condition knows which mesh nodes lie at its surface (checked by the assert),
// and it sets 0 in binaryMask at the appropriate location for these nodes, if the boundary
// condition type is currently Dirichlet.
void abstractBoundaryCondition::setDirichletConditionsIfNecessary(int* const binaryMask)
{
  assert(hasListOfMeshNodesAtThisBoundary);
  // set zero in the binaryMask at the locations necessary to impose Dirichlet at this surface
  for (auto node=listOfMeshNodesAtThisBoundary.begin(); node!=listOfMeshNodesAtThisBoundary.end(); node++)
  {
    binaryMask[*node] = 0;
  }
}