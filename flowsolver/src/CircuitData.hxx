#ifndef CIRCUITDATA_HXX_
#define CIRCUITDATA_HXX_

#include <map>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <cfloat>
#include <boost/shared_ptr.hpp>
#include "datatypesInCpp.hxx"
#include "gtest/gtest_prod.h"
#include "debuggingToolsForCpp.hxx"
#include "CircuitComponent.hxx"
#include "CircuitPressureNode.hxx"

class CircuitData
{
	friend class testMultidom;
	FRIEND_TEST(testMultidom,checkClosedDiodeWithRemainingOpenPathDetected);
	FRIEND_TEST(testMultidom,checkClosedDiodeWithoutRemainingOpenPathDetected);
public:
	CircuitData(const int hstep)
	: m_hstep(hstep)
	{
		m_flowPermittedAcross3DInterface = true;
		m_boundaryConditionTypeHasJustChanged = false;
	}

	virtual ~CircuitData()
	{
	}
	std::vector<boost::shared_ptr<CircuitComponent>> components;
	int index;

	int getLengthOfHistoryData()
	{
		try {
			return components.at(0)->m_entireFlowHistory.size();
		} catch (const std::exception& e) {
		    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
		    throw;
		}
	}
	
	// Begin metadata, updated with rebuildCircuitMetadata.
	int numberOfPrescribedPressures;
	int numberOfPrescribedFlows;
	int numberOfPressureNodes;
	int numberOfComponents;
	int m_numberOfVolumeTrackingComponenets;

	// These maps have indexInInputData as the key, and a shared_ptr to the relevant node/circuit as the mapped value.
	// Although this provides useful random access by indexInInputData, it is often useful to just use the std::map iterator to process them all.
	std::map<int,boost::shared_ptr<CircuitPressureNode>> mapOfPressureNodes; // Utility data structure, containing all the pressure nodes of the circuit, exactly once each.
	std::map<int,boost::shared_ptr<CircuitPressureNode>> mapOfPrescribedPressureNodes;
	std::map<int,boost::shared_ptr<CircuitComponent>> mapOfComponents; // Utility data structure containing all the component indices of the circuit, exactly once each.
	std::map<int,boost::shared_ptr<CircuitComponent>> mapOfPrescribedFlowComponents;
	std::map<int,boost::shared_ptr<CircuitComponent>> mapOfVolumeTrackingComponents;
	std::map<int,boost::shared_ptr<CircuitComponent>> mapOfPrescribedVolumeTrackingComponents;
	// End of medatata
	
	void rebuildCircuitMetadata();
	bool connectsTo3DDomain() const;
	virtual void initialiseNodeAndComponentAtInterface(int threeDInterfaceNodeIndex);
	virtual bool hasPrescribedFlowAcrossInterface() const;
	virtual bool hasPrescribedPressureAcrossInterface() const;
	void setupComponentNeighbourPointers();
	void switchDiodeStatesIfNecessary();
	void detectWhetherClosedDiodesStopAllFlowAt3DInterface();
	bool flowPermittedAcross3DInterface() const;
	bool boundaryConditionTypeHasJustChanged();

	std::vector<std::pair<int,double*>> getComponentInputDataIndicesAndFlows() const;
	std::vector<std::pair<int,double*>> getNodeInputDataIndicesAndPressures() const;
	std::vector<std::pair<int,double*>> getVolumeTrackingComponentInputDataIndicesAndVolumes() const;

	void setIndexOfNodeAtInterface(std::vector<int> indexToSet);
	int getIndexOfNodeAtInterface();
	int getIndexOfComponentConnectingToNodeAtInterface();

	void closeAllDiodes();
	boost::shared_ptr<CircuitComponent> getComponentByInputDataIndex(const int componentIndex);
	boost::shared_ptr<CircuitPressureNode> getNodeByInputDataIndex(const int componentIndex);

	boost::shared_ptr<CircuitPressureNode> ifExistsGetNodeOtherwiseConstructNode(const int indexInInputData_in, const circuit_nodal_pressure_prescription_t typeOfPrescribedPressure, const boost::shared_ptr<CircuitComponent> componentNeighbouringThisNode);
	std::vector<double*> getCapacitorNodalHistoryPressurePointers() const;
	double getSignForPrescribed3DInterfaceFlow() const;
	boost::shared_ptr<std::vector<std::pair<parameter_controller_t, int>>> getControlTypesAndComponentIndices() const;
protected:
	bool m_flowPermittedAcross3DInterface;
	std::vector<int> m_indexOfNodeAt3DInterface;
	void setIndicesOfNodesAtInterface(std::vector<int> indicesToSet);
	int toOneIndexing(const int oneIndexedValue);
	double m_signForPrescribed3DInterfaceFlow; // Necessary for if this component is at the 3D interface. If it's been connected to the 3D domain by the interfacing component's end-node, we need to switch the sign of the flow before prescribing it in the linear system for this boundary.
private:
	void rebuildCircuitPressureNodeMap();
	void switchBetweenDirichletAndNeumannCircuitDesign();
	int m_hstep;
	bool m_boundaryConditionTypeHasJustChanged;
};

class Netlist3DDomainReplacementCircuitData : public CircuitData
{
public:
	Netlist3DDomainReplacementCircuitData(const int hstep, const int numberOfNetlistsUsedAsBoundaryConditions)
	: CircuitData(hstep),
	m_numberOfNetlistsUsedAsBoundaryConditions(numberOfNetlistsUsedAsBoundaryConditions)
	{
	}
	bool hasPrescribedFlowAcrossInterface() const;
	bool hasPrescribedPressureAcrossInterface() const;
	void initialiseNodesAndComponentsAtInterface_vector(std::vector<int> threeDInterfaceNodeIndices);
	void setBoundaryPrescriptionsAndBoundaryConditionTypes(std::vector<std::pair<boundary_data_t,double>>& boundaryFlowsOrPressuresAsAppropriate);
	boost::shared_ptr<CircuitComponent> getDpDqResistorByIndex(int index);
	void addToMapOfDpDqResistors(int indexOfResistor, boost::shared_ptr<CircuitComponent> dpDqResistor);
	bool isADpDqResistor(const int componentIndex);
private:
	bool isNodeAtBoundaryInterface(int nodeIndex);
	const int m_numberOfNetlistsUsedAsBoundaryConditions;
	void givePrescribedPressureToBoundaryNode(int nodeIndex, double prescribedPressure);
	void givePrescribedFlowToBoundaryComponent(int componentIndex, double prescribedFlow);
	std::map<int,boost::shared_ptr<CircuitComponent>> m_mapOfDpDqResistors;

};

class ClosedLoopDownstreamCircuitData : public CircuitData
{
public:
	ClosedLoopDownstreamCircuitData(const int hstep) :
	CircuitData(hstep)
	{
	}
	void initialiseNodeAndComponentAtInterface(int threeDInterfaceNodeIndex)
	{
		// Does nothing, as there's no node at 3D interface. \todo this could be tidier, but not essential.
	}

};

#endif
