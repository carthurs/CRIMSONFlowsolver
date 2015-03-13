#ifndef CIRCUITDATA_HXX_
#define CIRCUITDATA_HXX_

#include <map>
#include <vector>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include "datatypesInCpp.hxx"
#include "gtest/gtest_prod.h"
#include "debuggingToolsForCpp.hxx"

// Forward declaration:
class CircuitPressureNode;

class CircuitComponent
{
public:
	boost::shared_ptr<CircuitPressureNode> startNode;
	// bool m_startNodeConnectsToDiode;
	// boost::shared_ptr<CircuitPressureNode> startNodeIfNeighbouringDiodeExistsAndIsOpen;

	boost::shared_ptr<CircuitPressureNode> endNode;
	// bool m_endNodeConnectsToDiode;
	// boost::shared_ptr<CircuitPressureNode> endNodeIfNeighbouringDiodeExistsAndIsOpen;

	std::vector<boost::weak_ptr<CircuitComponent>> neighbouringComponentsAtStartNode;
	std::vector<boost::weak_ptr<CircuitComponent>> neighbouringComponentsAtEndNode;

	int prescribedFlowPointerIndex;
	
	double parameterValueFromInputData; // for diodes only. Stores a value from netlist_surfaces.dat to be set as the currentParameterValue (resistance) when the diode is open.
	int indexInInputData;
	int indexLocalToSubcircuit;
	circuit_component_flow_prescription_t prescribedFlowType;
	double valueOfPrescribedFlow;
	double m_signForPrescribed3DInterfaceFlow; // Necessary for if this component is at the 3D interface. If it's been connected to the 3D domain by its end-node, we need to switch the sign of the flow before prescribing it in the linear system for this boundary.
	double flow;
	double historyFlow;
	bool hasHistoryFlow;
	bool hasTrackedVolume;
	bool hasHistoryVolume;
	bool permitsFlow; // for diodes in particular
	std::vector<double> m_entireFlowHistory;
	CircuitComponent(const int hstep, const bool thisIsARestartedSimulation)
	: m_hstep(hstep),
	  m_thisIsARestartedSimulation(thisIsARestartedSimulation)
	{
		m_type = Component_Null;
		prescribedFlowType = Flow_Null;
		hasHistoryFlow = false;
		m_entireFlowHistory.reserve(m_hstep);
		if (m_thisIsARestartedSimulation)
		{
			bool fixThisForRestart=false;
            assert(fixThisForRestart);
            flow = -1.0;
            permitsFlow = true;
            m_connectsToNodeAtInterface = false;
            hasTrackedVolume = false;
			hasHistoryVolume = false;
		}
		else
		{
			flow = 0.0;
			permitsFlow = true;
			m_connectsToNodeAtInterface = false;
			hasTrackedVolume = false;
			hasHistoryVolume = false;
		}
		prescribedFlowPointerIndex = 0;
	}

	virtual ~CircuitComponent()
	{
	}

	bool hasNonnegativePressureGradientOrForwardFlow(); // whether the diode should be open
	bool connectsToNodeAtInterface();
	void setConnectsToNodeAtInterface();
	void enableDiodeFlow();
	void disableDiodeFlow();
	circuit_component_t& getType();
	double* getParameterPointer();
	void setParameterValue(double const parameterValue);
private:
	circuit_component_t m_type;
	double m_currentParameterValue; // resistance or compliance or inductance or elastance etc.

	const int m_hstep;
	const bool m_thisIsARestartedSimulation;
	bool m_connectsToNodeAtInterface;
};

// A slightly more complicated class of component, which prescribes its pressure
// in the circuit depending on its elastance and its stored volume.
// Think of it more as a chamber than as a node.
class VolumeTrackingPressureChamber : public CircuitComponent
{
public:
	VolumeTrackingPressureChamber(const int hstep, const bool thisIsARestartedSimulation)
	: CircuitComponent(hstep, thisIsARestartedSimulation)
	{
		assert(!thisIsARestartedSimulation);
		m_storedVolume = 130000.0; // default; can be changed later if necessary
		m_unstressedVolume = 0.0; // default; can be changed later if necessary
		m_entireVolumeHistory.reserve(hstep);
		m_enforceZeroVolumePrescription = false;
	}

	void recordVolumeInHistory()
	{
		m_entireVolumeHistory.push_back(m_storedVolume);
	}

	double getVolumeHistoryAtTimestep(int timestep)
	{
		return m_entireVolumeHistory.at(timestep);
	}

	void setStoredVolume(const double newVolume)
	{
		m_storedVolume = newVolume;
	}
	// The /proposed/ volume is the one which gets checked for negative (invalid) values
	// so that we can detect such invalid cases, and take steps to remedy.
	void setProposedVolume(const double proposedVolume)
	{
		m_proposedVolume = proposedVolume;
	}
	double getVolume()
	{
		return m_storedVolume;
	}
	double getProposedVolume()
	{
		std::cout<<"proposed volume was: " << m_proposedVolume << std::endl;
		return m_proposedVolume;	
	}
	double getHistoryVolume()
	{
		return m_historyVolume;
	}
	void cycleHistoryVolume()
	{
		m_historyVolume = m_storedVolume;
	}
	void passPressureToStartNode();

	double* getPointerToElastance()
	{
		return &m_currentParameterValue;
	}

	double getElastance()
	{
		return m_currentParameterValue;
	}

	bool zeroVolumeShouldBePrescribed()
	{
		return m_enforceZeroVolumePrescription;
	}

	void enforceZeroVolumePrescription()
	{
		m_enforceZeroVolumePrescription = true;
	}

	void resetZeroVolumePrescription()
	{
		m_enforceZeroVolumePrescription = false;
	}
private:
	double m_pressure;
	double m_storedVolume;
	double m_proposedVolume; // this holds volumes temporarily, so that we can check them for invalid negative values & do something about it if necessary
	double m_unstressedVolume;
	double m_historyVolume;
	bool m_enforceZeroVolumePrescription;
	std::vector<double> m_entireVolumeHistory;
};

class CircuitPressureNode
{
public:
	double historyPressure;
	bool hasHistoryPressure;
	circuit_nodal_pressure_prescription_t prescribedPressureType;
	int prescribedPressurePointerIndex;
	int indexInInputData;
	int indexLocalToSubcircuit;
	std::vector<double> m_entirePressureHistory;

	std::vector<boost::weak_ptr<CircuitComponent>> listOfComponentstAttachedToThisNode;
	CircuitPressureNode(const int indexInInputData_in, const circuit_nodal_pressure_prescription_t typeOfPrescribedPressure, const int hstep)
	: indexInInputData(indexInInputData_in),
	prescribedPressureType(typeOfPrescribedPressure),
	m_hstep(hstep)
	{
		hasHistoryPressure = false;
	    m_isAtBoundary = false;
	    m_entirePressureHistory.reserve(m_hstep);
	    prescribedPressurePointerIndex = 0;
	}

	void setIsAtBoundary()
	{
		m_isAtBoundary = true;
	}

	bool isAtBoundary()
	{
		return m_isAtBoundary;
	}

	double getPressure()
	{
		return pressure;
	}
	void setPressure(const double pressure_in)
	{
		pressure = pressure_in;
	}
protected:
	double pressure;
	const int m_hstep;
private:
	bool m_isAtBoundary;
};

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
		return components.at(0)->m_entireFlowHistory.size();
	}

	// // copy constructor
	// CircuitData(const CircuitData &sourceCircuitData)
	// : m_hstep(sourceCircuitData.m_hstep),
	// components(sourceCircuitData.components),
	// index(sourceCircuitData.index),
	// numberOfPrescribedPressures(sourceCircuitData.numberOfPrescribedPressures),
	// numberOfPrescribedFlows(sourceCircuitData.numberOfPrescribedFlows),
	// numberOfPressureNodes(sourceCircuitData.numberOfPressureNodes),
	// numberOfComponents(sourceCircuitData.numberOfComponents),
	// mapOfPressureNodes(sourceCircuitData.mapOfPressureNodes),
	// mapOfPrescribedPressureNodes(sourceCircuitData.mapOfPrescribedPressureNodes),
	// mapOfComponents(sourceCircuitData.mapOfComponents),
	// mapOfPrescribedFlowComponents(sourceCircuitData.mapOfPrescribedFlowComponents)
	// {
	// }

	// // assignment operator
	// CircuitData operator=(CircuitData rhs)
	// {
	// 	return CircuitData(rhs);
	// }
	
	// Begin metadata, updated with rebuildCircuitMetadata.
	int numberOfPrescribedPressures;
	int numberOfPrescribedFlows;
	int numberOfPressureNodes;
	int numberOfComponents;
	int m_numberOfVolumeTrackingPressureChambers;

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
	void generateNodeAndComponentIndicesLocalToSubcircuit();
	virtual void initialiseNodeAndComponentAtInterface(int threeDInterfaceNodeIndex);
	virtual bool hasPrescribedFlowAcrossInterface() const;
	virtual bool hasPrescribedPressureAcrossInterface() const;
	void setupComponentNeighbourPointers();
	void switchDiodeStatesIfNecessary();
	void detectWhetherClosedDiodesStopAllFlowAt3DInterface();
	bool flowPermittedAcross3DInterface() const;
	bool boundaryConditionTypeHasJustChanged();

	void setIndexOfNodeAtInterface(std::vector<int> indexToSet);
	int getIndexOfNodeAtInterface();
	int getIndexOfComponentConnectingToNodeAtInterface();

	void closeAllDiodes();
	boost::shared_ptr<CircuitComponent> getComponentByInputDataIndex(const int componentIndex);

	boost::shared_ptr<CircuitPressureNode> ifExistsGetNodeOtherwiseConstructNode(const int indexInInputData_in, const circuit_nodal_pressure_prescription_t typeOfPrescribedPressure, const boost::shared_ptr<CircuitComponent> componentNeighbouringThisNode);

protected:
	bool m_flowPermittedAcross3DInterface;
	std::vector<int> m_indexOfNodeAt3DInterface;
	void setIndicesOfNodesAtInterface(std::vector<int> indicesToSet);
	int toOneIndexing(const int oneIndexedValue);
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
private:
	bool isNodeAtBoundaryInterface(int nodeIndex);
	const int m_numberOfNetlistsUsedAsBoundaryConditions;
	void givePrescribedPressureToBoundaryNode(int nodeIndex, double prescribedPressure);
	void givePrescribedFlowToBoundaryComponent(int componentIndex, double prescribedFlow);

};

#endif
