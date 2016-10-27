#ifndef CIRCUITCOMPONENT_HXX_
#define CIRCUITCOMPONENT_HXX_

#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "datatypesInCpp.hxx"

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

	boost::shared_ptr<CircuitPressureNode> getStartNode();
	boost::shared_ptr<CircuitPressureNode> getEndNode();

	std::vector<boost::weak_ptr<CircuitComponent>> neighbouringComponentsAtStartNode;
	std::vector<boost::weak_ptr<CircuitComponent>> neighbouringComponentsAtEndNode;

	int prescribedFlowPointerIndex;
	
	double parameterValueFromInputData; // for diodes only. Stores a value from netlist_surfaces.dat to be set as the currentParameterValue (resistance) when the diode is open.
	circuit_component_flow_prescription_t prescribedFlowType;
	double flow;
	double historyFlow;
	bool hasHistoryFlow;
	// bool hasTrackedVolume;
	CircuitComponent(const int hstep, const bool thisIsARestartedSimulation)
	: m_hstep(hstep),
	  m_thisIsARestartedSimulation(thisIsARestartedSimulation)
	{
		m_type = Component_Null;
		prescribedFlowType = Flow_Null;
		hasHistoryFlow = false;
		m_entireFlowHistory.reserve(m_hstep);
		// if (m_thisIsARestartedSimulation)
		// {
		// 	bool fixThisForRestart=false;
  //           assert(fixThisForRestart);
  //           flow = -1.0;
  //           m_permitsFlow = true;
  //           m_connectsToNodeAtInterface = false;
  //           m_hasPythonParameterController = false;
  //           // hasTrackedVolume = false;
		// 	// hasHistoryVolume = false;
		// 	m_hasPrescribedFlow = false;
		// }
		// else
		// {
			flow = 0.0;
			m_permitsFlow = true;
			m_connectsToNodeAtInterface = false;
			m_hasPythonParameterController = false; //\todo fix this for restart
			// hasTrackedVolume = false;
			m_hasHistoryVolume = false;
			m_hasPrescribedFlow = false;
		// }
		prescribedFlowPointerIndex = 0;

		// m_flowHistoryBufferSize = 100;
	    // m_flowHistoryBuffer.resize(m_flowHistoryBufferSize);
		// m_flowHistoryBufferNextWriteIndex = 0;
	}

	virtual ~CircuitComponent()
	{
	}

	bool hasUserDefinedExternalPythonScriptParameterController() const;
	std::string getPythonControllerName(const parameter_controller_t controllerType) const;
	void addPythonControllerName(const parameter_controller_t controlType, const std::string pythonParameterControllerName);

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
	double getPrescribedFlow() const;
	void setPrescribedFlow(const double prescribedFlow);
	void setHasNoPrescribedFlow();
	double* getPointerToFixedFlowPrescription();
	bool hasPrescribedFlow() const;
	void setRestartFlowFromHistory();
	void setHasHistoryVolume(const bool hasHistoryVolume);
	bool getHasHistoryVolume();
	void appendToFlowHistory(const double flow);
	int getFlowHistoryLength() const;
	double getFromFlowHistoryByTimestepIndex(const int timestepIndex) const;
	parameter_controller_t getControlType() const;
protected:
	double m_currentParameterValue; // resistance or compliance or inductance or elastance etc.
	bool m_hasPythonParameterController;
	std::map<parameter_controller_t, std::string> m_pythonParameterControllerNames;
private:
	circuit_component_t m_type;
	bool m_permitsFlow; // for diodes in particular
	double m_valueOfPrescribedFlow;
	bool m_hasHistoryVolume;

	std::vector<double> m_entireFlowHistory;
	
	int m_indexInInputData;
	const int m_hstep;
	const bool m_thisIsARestartedSimulation;
	bool m_connectsToNodeAtInterface;

	bool m_hasPrescribedFlow;

	// std::vector<double> m_flowHistoryBuffer;
	// int m_flowHistoryBufferNextWriteIndex;
	// int m_flowHistoryBufferSize;
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
		//assert(!thisIsARestartedSimulation);
		m_entireVolumeHistory.reserve(hstep);
		m_enforceZeroVolumePrescription = false;
	}

	double getVolumeHistoryAtTimestep(int timestep);
	void setVolumeHistoryAtTimestep(double historyVolume);
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
	void setRestartVolumeFromHistory();
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
	VolumeTrackingPressureChamber(const int hstep, const bool thisIsARestartedSimulation, const double initialVolume, const double initialUnstressedVolume)
	: VolumeTrackingComponent(hstep, thisIsARestartedSimulation, initialVolume)
	{
		//assert(!thisIsARestartedSimulation);
		m_unstressedVolume = initialUnstressedVolume;
	}
	
	void setStoredVolume(const double newVolume);
	double getUnstressedVolume();
	double* getUnstressedVolumePointer();
private:
	double m_unstressedVolume;
	void passPressureToStartNode();
};

#endif