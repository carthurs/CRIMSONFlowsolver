#ifndef NETLISTCLOSEDLOOPDOWNSTREAMCIRCUIT_HXX_
#define NETLISTCLOSEDLOOPDOWNSTREAMCIRCUIT_HXX_

#include "NetlistCircuit.hxx"

class NetlistClosedLoopDownstreamCircuit : public NetlistCircuit
{
public:
	NetlistClosedLoopDownstreamCircuit(const int hstep, const bool thisIsARestartedSimulation, const double alfi, const double delt)
	: NetlistCircuit(hstep, thisIsARestartedSimulation, alfi, delt)
	{
		m_downstreamCircuitIndex = s_numberOfDownstreamCircuits;
		s_numberOfDownstreamCircuits++;

		std::stringstream pressureFileNameBuilder;
		pressureFileNameBuilder << "netlistPressures_downstreamCircuit_" << m_downstreamCircuitIndex << ".dat";
		m_PressureHistoryFileName = pressureFileNameBuilder.str();

		std::stringstream flowFileNameBuilder;
		flowFileNameBuilder << "netlistFlows_downstreamCircuit_" << m_downstreamCircuitIndex << ".dat";
		m_FlowHistoryFileName = flowFileNameBuilder.str();

		std::stringstream volumeFileNameBuilder;
		volumeFileNameBuilder << "netlistVolumes_downstreamCircuit_" << m_downstreamCircuitIndex << ".dat";
		m_VolumeHistoryFileName = volumeFileNameBuilder.str();
	}

	void createCircuitDescription();
	void initialiseCircuit();

	void getMatrixContribution(const double alfi_delt, Mat& matrixFromThisDownstreamCircuit);
	void getRHSContribution(const int timestepNumber, Vec& rhsFromThisDownstreamCircuit);

	int getCircuitIndex() const;

	int convertInterfaceNodeIndexFromDownstreamToUpstreamCircuit(const int sharedNodeDownstreamIndex) const;

	void getSharedNodeDownstreamAndUpstreamAndCircuitUpstreamIndices(std::vector<int>& downstreamNodeIndices, std::vector<int>& upstreamNodeIndices, std::vector<int>& upstreamSurfaceIndices) const;

	~NetlistClosedLoopDownstreamCircuit()
	{
		s_numberOfDownstreamCircuits--;
	}
private:
	static int s_numberOfDownstreamCircuits;
	int m_downstreamCircuitIndex;
	int m_numberOfConnectedBoundaryConditions;
	
	std::vector<int> m_connectedCircuitSurfaceIndices;
	std::vector<int> m_localInterfacingNodes;
	std::vector<int> m_remoteInterfacingNodes;

	std::set<int> m_pressureNodesWhichConnectToBoundaryCircuits;
	std::map<int,int> m_circuitInterfaceNodeIndexMapDownstreamToUpstream;

	int m_numberOfNodesConnectingToAnotherCircuit;

	void appendClosedLoopSpecificCircuitDescription();
	bool kirchoffEquationAtNodeDeferredToInterfacingCircuit(const int nodeIndex) const;

	// Disabling methods that should never be called:
	void detectWhetherClosedDiodesStopAllFlowAt3DInterface();
};

#endif