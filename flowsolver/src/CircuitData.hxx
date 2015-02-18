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
	circuit_component_t type;
	boost::shared_ptr<CircuitPressureNode> startNode;
	// bool m_startNodeConnectsToDiode;
	// boost::shared_ptr<CircuitPressureNode> startNodeIfNeighbouringDiodeExistsAndIsOpen;

	boost::shared_ptr<CircuitPressureNode> endNode;
	// bool m_endNodeConnectsToDiode;
	// boost::shared_ptr<CircuitPressureNode> endNodeIfNeighbouringDiodeExistsAndIsOpen;

	std::vector<boost::weak_ptr<CircuitComponent>> neighbouringComponentsAtStartNode;
	std::vector<boost::weak_ptr<CircuitComponent>> neighbouringComponentsAtEndNode;
	
	double currentParameterValue; // resistance or compliance or inductance etc.
	double parameterValueFromInputData; // for diodes only. Stores a value from netlist_surfaces.dat to be set as the currentParameterValue (resistance) when the diode is open.
	int indexInInputData;
	int indexLocalToSubcircuit;
	circuit_component_flow_prescription_t prescribedFlowType;
	double valueOfPrescribedFlow;
	double signForPrescribed3DInterfaceFlow; // Necessary for if this component is at the 3D interface. If it's been connected to the 3D domain by its end-node, we need to switch the sign of the flow before prescribing it in the linear system for this boundary.
	double flow;
	double historyFlow;
	bool hasHistoryFlow;
	bool permitsFlow; // for diodes in particular
	std::vector<double> m_entireFlowHistory;
	CircuitComponent(const int hstep, const bool thisIsARestartedSimulation)
	: m_hstep(hstep),
	  m_thisIsARestartedSimulation(thisIsARestartedSimulation)
	{
		type = Component_Null;
		prescribedFlowType = Flow_Null;
		hasHistoryFlow = false;
		m_entireFlowHistory.reserve(m_hstep);
		if (m_thisIsARestartedSimulation)
		{
			bool fixThisForRestart=false;
            assert(fixThisForRestart);
            flow = -1.0;
            permitsFlow = true;
            m_connectsToNodeAt3DInterface = false;
		}
		else
		{
			flow = 0.0;
			permitsFlow = true;
			m_connectsToNodeAt3DInterface = false;
		}
	}

	bool hasNonnegativePressureGradientOrForwardFlow(); // whether the diode should be open
	bool connectsToNodeAt3DInterface();
	void setConnectsToNodeAt3DInterface();
private:
	const int m_hstep;
	const bool m_thisIsARestartedSimulation;
	bool m_connectsToNodeAt3DInterface;
};

class CircuitPressureNode
{
public:
	double historyPressure;
	bool hasHistoryPressure;
	circuit_nodal_pressure_prescription_t prescribedPressureType;
	int indexInInputData;
	int indexLocalToSubcircuit;
	bool m_connectsTo3DDomain;
	std::vector<double> m_entirePressureHistory;

	std::vector<boost::weak_ptr<CircuitComponent>> listOfComponentstAttachedToThisNode;
	CircuitPressureNode(const int indexInInputData_in, const circuit_nodal_pressure_prescription_t typeOfPrescribedPressure, const int hstep)
	: indexInInputData(indexInInputData_in),
	prescribedPressureType(typeOfPrescribedPressure),
	m_hstep(hstep)
	{
		hasHistoryPressure = false;
	    m_connectsTo3DDomain = false;
	    m_entirePressureHistory.reserve(m_hstep);
	}
	virtual ~CircuitPressureNode()
	{
	}

	virtual double getPressure()
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
};

// A slightly more complicated class of node, which prescribes its pressure
// in the circuit depending on its compliance and its stored volume.
// Think of it more as a chamber than as a node.
class VolumeTrackingPressureChamber : public CircuitPressureNode
{
public:
	VolumeTrackingPressureChamber(const int indexInInputData_in, const circuit_nodal_pressure_prescription_t typeOfPrescribedPressure, const int hstep)
	: CircuitPressureNode(indexInInputData_in, typeOfPrescribedPressure, hstep)
	{
		prescribedPressureType = Pressure_VolumeDependent;
		storedVolume = 0.0; // default; can be changed later if necessary
		unstressedVolume = 0.0; // default; can be changed later if necessary
		compliance = 1.0;
		m_entireVolumeHistory.reserve(m_hstep);
	}

	std::vector<double> m_entireVolumeHistory;

	void updateStoredVolume(const double delt)
	{
		double volumeChange = 0.0;

		// Find the attached components and add their flow contributions to the volume in the VolumeTrackingPressureChamber:
		for (auto attachedComponent=listOfComponentstAttachedToThisNode.begin(); attachedComponent!=listOfComponentstAttachedToThisNode.end(); attachedComponent++)
		{
			// If the VolumeTrackingPressureChamber is the startNode of the current component, the flow
			// is /out/ of the chamber:
			if ((*attachedComponent).lock()->startNode->indexInInputData == indexInInputData)
			{
				volumeChange -= (*attachedComponent).lock()->flow * delt;
				std::cout << "just did start node flow " << (*attachedComponent).lock()->flow * delt << std::endl;
			}
			else if ((*attachedComponent).lock()->endNode->indexInInputData == indexInInputData) // If VolumeTrackingPressureChamber is the current component's endNode, then the flow is into the chamber.
			{
				volumeChange += (*attachedComponent).lock()->flow * delt;
				std::cout << "just did end node flow " << (*attachedComponent).lock()->flow * delt << std::endl;
			}
			else
			{
				throw std::logic_error("EE: Reached an impossible internal location. Please contact the developers.");
			}
		}

		storedVolume = storedVolume + volumeChange;
	}
	double getVolume()
	{
		return storedVolume;
	}
	double getPressure()
	{
		pressure = (storedVolume - unstressedVolume)/compliance;
		return pressure;
	}
private:
	double storedVolume;
	double compliance;
	double unstressedVolume;
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
	std::vector<boost::shared_ptr<CircuitComponent>> components;
	int index;

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
	// End of medatata
	
	void rebuildCircuitMetadata();
	bool connectsTo3DDomain() const;
	void generateNodeAndComponentIndicesLocalToSubcircuit();
	void initialiseNodeAndComponentAt3DInterface(int threeDInterfaceNodeIndex);
	void setupComponentNeighbourPointers();
	void switchDiodeStatesIfNecessary();
	void detectWhetherClosedDiodesStopAllFlowAt3DInterface();
	bool flowPermittedAcross3DInterface() const;
	bool boundaryConditionTypeHasJustChanged();

	void setIndexOfNodeAt3DInterface(int indexToSet);
	int getIndexOfNodeAt3DInterface();

	boost::shared_ptr<CircuitPressureNode> ifExistsGetNodeOtherwiseConstructNode(const int indexInInputData_in, const circuit_nodal_pressure_prescription_t typeOfPrescribedPressure, const boost::shared_ptr<CircuitComponent> componentNeighbouringThisNode);
private:
	int toOneIndexing(const int oneIndexedValue);
	void rebuildCircuitPressureNodeMap();
	void switchBetweenDirichletAndNeumannCircuitDesign();
	int m_hstep;
	bool m_flowPermittedAcross3DInterface;
	bool m_boundaryConditionTypeHasJustChanged;
	int m_indexOfNodeAt3DInterface;
};

#endif
