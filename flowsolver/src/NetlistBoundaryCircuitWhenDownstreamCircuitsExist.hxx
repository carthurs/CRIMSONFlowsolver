#ifndef NETLISTBOUNDARYCIRCUITWHENDOWNSTREAMCIRCUITSEXIST_HXX_
#define NETLISTBOUNDARYCIRCUITWHENDOWNSTREAMCIRCUITSEXIST_HXX_

#include "NetlistCircuit.hxx"

// Forward declarations:
class ClosedLoopDownstreamSubsection;

class NetlistBoundaryCircuitWhenDownstreamCircuitsExist : public NetlistCircuit
{
public:
	NetlistBoundaryCircuitWhenDownstreamCircuitsExist(const int hstep, const int surfaceIndex, const int indexOfThisNetlistLPN, const bool thisIsARestartedSimulation, const double alfi, const double delt, const std::vector<boost::weak_ptr<ClosedLoopDownstreamSubsection>> downstreamSubcircuits, const int startingTimestepIndex)
	:NetlistCircuit(hstep, surfaceIndex, indexOfThisNetlistLPN, thisIsARestartedSimulation, alfi, delt, startingTimestepIndex),
	m_netlistDownstreamLoopClosingSubcircuits(downstreamSubcircuits)
	{}
	void initialiseCircuit();
	std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timeAtStepNplus1, const double alfi_delt);
	void getMatrixContribution(const double alfi_delt, Mat& matrixFromThisBoundary);
	void getRHSContribution(const int timestepNumber, Vec& rhsFromThisBoundary);
	int getLocationOf3DInterfaceComputedFlowInSolutionVector() const;
	int getLocationOf3DInterfaceComputedPressureInSolutionVector() const;
	int getColumnOf3DInterfacePrescribedPressureInLinearSystem() const;
	int getColumnOf3DInterfacePrescribedFlowInLinearSystem() const;
	std::pair<boundary_data_t,double> computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement(const int timestepNumber);
	int getCircuitIndex() const;
	void updateLPN(const int timestepNumber);
protected:
private:
	std::set<int> m_pressureNodesWhichConnectToDownstreamCircuits;
	int m_numberOfNodesConnectingToAnotherCircuit;
	std::vector<boost::weak_ptr<ClosedLoopDownstreamSubsection>> m_netlistDownstreamLoopClosingSubcircuits;
	bool kirchoffEquationAtNodeDeferredToInterfacingCircuit(const int nodeIndex) const;
};

#endif