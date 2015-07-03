#ifndef CIRCUITDATA_HXX_
#define CIRCUITDATA_HXX_

#include <map>
#include <unordered_map>
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
	circuit_component_flow_prescription_t prescribedFlowType;
	double valueOfPrescribedFlow;
	double m_signForPrescribed3DInterfaceFlow; // Necessary for if this component is at the 3D interface. If it's been connected to the 3D domain by its end-node, we need to switch the sign of the flow before prescribing it in the linear system for this boundary.
	double flow;
	double historyFlow;
	bool hasHistoryFlow;
	bool hasTrackedVolume;
	bool hasHistoryVolume;
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
            m_permitsFlow = true;
            m_connectsToNodeAtInterface = false;
            m_hasPythonParameterController = false;
            hasTrackedVolume = false;
			hasHistoryVolume = false;
		}
		else
		{
			flow = 0.0;
			m_permitsFlow = true;
			m_connectsToNodeAtInterface = false;
			m_hasPythonParameterController = false; //\todo fix this for restart
			hasTrackedVolume = false;
			hasHistoryVolume = false;
		}
		prescribedFlowPointerIndex = 0;
	}

	virtual ~CircuitComponent()
	{
	}

	bool hasUserDefinedExternalPythonScriptParameterController() const;
	std::string getPythonControllerName() const;
	void setPythonControllerName(const std::string pythonParameterControllerName);

	bool hasNonnegativePressureGradientOrForwardFlow(); // whether the diode should be open
	bool connectsToNodeAtInterface();
	void setConnectsToNodeAtInterface();
	void enableDiodeFlow();
	void disableDiodeFlow();
	circuit_component_t& getType();
	double* getParameterPointer();
	double* getFlowPointer();
	void setParameterValue(double const parameterValue);
	int getIndex() const;
	void setIndex(const int index);
	bool permitsFlow() const;
protected:
	double m_currentParameterValue; // resistance or compliance or inductance or elastance etc.
	bool m_hasPythonParameterController;
	std::string m_pythonParameterControllerName;
private:
	circuit_component_t m_type;
	bool m_permitsFlow; // for diodes in particular
	
	int m_indexInInputData;
	const int m_hstep;
	const bool m_thisIsARestartedSimulation;
	bool m_connectsToNodeAtInterface;
};

// A slightly more complicated class of component, which prescribes its pressure
// in the circuit depending on its elastance and its stored volume.
// Think of it more as a chamber than as a node.
class VolumeTrackingComponent : public CircuitComponent
{
public:
	VolumeTrackingComponent(const int hstep, const bool thisIsARestartedSimulation, const double initialVolume)
	: CircuitComponent(hstep, thisIsARestartedSimulation),
	m_storedVolume(initialVolume)
	{
		assert(!thisIsARestartedSimulation);
		m_entireVolumeHistory.reserve(hstep);
		m_enforceZeroVolumePrescription = false;
	}

	double getVolumeHistoryAtTimestep(int timestep);
	virtual void setStoredVolume(const double newVolume);
	void setProposedVolume(const double proposedVolume);
	double getVolume();
	double* getVolumePointer();
	double getProposedVolume();
	double getHistoryVolume();
	void cycleHistoryVolume();
	double getElastance();
	bool zeroVolumeShouldBePrescribed();
	void enforceZeroVolumePrescription();
	void resetZeroVolumePrescription();
	void recordVolumeInHistory();
protected:
	double m_pressure;
	double m_storedVolume;
	double m_proposedVolume; // this holds volumes temporarily, so that we can check them for invalid negative values & do something about it if necessary
	double m_historyVolume;
	bool m_enforceZeroVolumePrescription;
	std::vector<double> m_entireVolumeHistory;
};

class VolumeTrackingPressureChamber : public VolumeTrackingComponent
{
public:
	VolumeTrackingPressureChamber(const int hstep, const bool thisIsARestartedSimulation, const double initialVolume)
	: VolumeTrackingComponent(hstep, thisIsARestartedSimulation, initialVolume)
	{
		assert(!thisIsARestartedSimulation);
		m_unstressedVolume = 0.0; // default; can be changed later if necessary
	}
	
	void setStoredVolume(const double newVolume);
private:
	double m_unstressedVolume;
	void passPressureToStartNode();
};

class CircuitPressureNode
{
public:
	double historyPressure;
	bool hasHistoryPressure;
	circuit_nodal_pressure_prescription_t prescribedPressureType;
	int prescribedPressurePointerIndex;
	std::vector<double> m_entirePressureHistory;

	std::vector<boost::weak_ptr<CircuitComponent>> listOfComponentstAttachedToThisNode;
	CircuitPressureNode(const int indexInInputData, const circuit_nodal_pressure_prescription_t typeOfPrescribedPressure, const int hstep)
	: m_indexInInputData(indexInInputData),
	prescribedPressureType(typeOfPrescribedPressure),
	m_hstep(hstep)
	{
		hasHistoryPressure = false;
	    m_isAtBoundary = false;
	    m_hasPythonParameterController = false;
	    m_entirePressureHistory.reserve(m_hstep);
	    prescribedPressurePointerIndex = 0;
	}

	double* getPressurePointer();
	double* getPointerToFixedPressurePrescription();
	void setIsAtBoundary();
	bool isAtBoundary() const;
	double getPressure();
	void setPressure(const double pressure_in);
	void setPrescribedPressure(const double prescribedPressure);
	int getIndex() const;
	void setPythonControllerName(const std::string pythonParameterControllerName);
	bool hasUserDefinedExternalPythonScriptParameterController() const;
	std::string getPythonControllerName() const;

protected:
	double pressure;
	const int m_hstep;
private:
	bool m_isAtBoundary;
	const int m_indexInInputData;
	bool m_hasPythonParameterController;
	std::string m_pythonParameterControllerName;
	double m_fixedPressure;
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
