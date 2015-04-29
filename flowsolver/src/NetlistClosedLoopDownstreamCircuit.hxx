#ifndef NETLISTCLOSEDLOOPDOWNSTREAMCIRCUIT_HXX_
#define NETLISTCLOSEDLOOPDOWNSTREAMCIRCUIT_HXX_

#include "NetlistCircuit.hxx"

class twoDimensionalMap
{
public:
	void insert(const int key1, const int key2, const int mappedValue)
	{
		std::map<int,int> toInsert;
		toInsert.insert(std::make_pair(key2, mappedValue));
		m_data.insert(std::make_pair(key1, toInsert));
	}

	int at(const int key1, const int key2) const
	{
		return m_data.at(key1).at(key2);
	}

	int size() const
	{
		int sizeToReturn = 0;
		for (auto firstKeyLevel = m_data.begin(); firstKeyLevel != m_data.end(); firstKeyLevel++)
		{
			sizeToReturn += firstKeyLevel->second.size();
		}
		return sizeToReturn;
	}

private:
	std::map<int,std::map<int,int>> m_data;
};

class NetlistClosedLoopDownstreamCircuit : public NetlistCircuit
{
public:
	NetlistClosedLoopDownstreamCircuit(const int hstep, const bool thisIsARestartedSimulation, const double alfi, const double delt)
	: NetlistCircuit(hstep, s_numberOfDownstreamCircuits, thisIsARestartedSimulation, alfi, delt)
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

		mp_circuitData = boost::shared_ptr<CircuitData> (new ClosedLoopDownstreamCircuitData(hstep));
	}

	void createCircuitDescription();
	void initialiseCircuit();

	void getMatrixContribution(const double alfi_delt, Mat& matrixFromThisDownstreamCircuit);
	void getRHSContribution(const int timestepNumber, Vec& rhsFromThisDownstreamCircuit);
	void initialiseAtStartOfTimestep();

	int getCircuitIndex() const;

	std::set<int> getPressureNodesConnectingToUpstreamCircuits() {return m_pressureNodesWhichConnectToBoundaryCircuits;};

	int convertInterfaceNodeIndexFromDownstreamToUpstreamCircuit(const int surfaceIndex, const int sharedNodeDownstreamIndex) const;
	int convertInterfaceNodeIndexFromUpstreamToDownstreamCircuit(const int sharedNodeUpstreamIndex) const;
	bool boundaryConditionCircuitConnectsToThisDownstreamSubsection(const int boundaryConditionIndex) const;

	void getSharedNodeDownstreamAndUpstreamAndCircuitUpstreamIndices(std::vector<int>& downstreamNodeIndices, std::vector<int>& upstreamNodeIndices, std::vector<int>& upstreamSurfaceIndices) const;
	void giveNodesAndComponentsTheirUpdatedValuesFromSolutionVector(const std::vector<PetscScalar> solutionEntriesForDownstreamCircuit);

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
	twoDimensionalMap m_circuitInterfaceNodeIndexMapDownstreamToUpstream;
	std::map<int,int> m_circuitInterfaceNodeIndexMapUpstreamToDownstream;

	std::set<int> m_setOfAttachedBoundaryConditionIndices;

	int m_numberOfNodesConnectingToAnotherCircuit;

	void appendClosedLoopSpecificCircuitDescription();
	bool kirchoffEquationAtNodeDeferredToInterfacingCircuit(const int nodeIndex) const;

	// Disabling methods that should never be called:
	void detectWhetherClosedDiodesStopAllFlowAt3DInterface();
};

#endif