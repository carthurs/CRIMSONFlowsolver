#ifndef CIRCUITPRESSURENODE_HXX_
#define CIRCUITPRESSURENODE_HXX_

#include <vector>
#include <string>

#include <boost/weak_ptr.hpp>

#include "datatypesInCpp.hxx"
#include "CircuitComponent.hxx"

class CircuitPressureNode
{
public:
	std::vector<boost::weak_ptr<CircuitComponent>> listOfComponentstAttachedToThisNode;
	CircuitPressureNode(const int indexInInputData, const circuit_nodal_pressure_prescription_t typeOfPrescribedPressure, const int hstep)
	: m_indexInInputData(indexInInputData),
	m_prescribedPressureType(typeOfPrescribedPressure),
	m_hstep(hstep)
	{
		m_hasHistoryPressure = false;
	    m_isAtBoundary = false;
	    m_hasPythonParameterController = false;
	    m_entirePressureHistory.reserve(m_hstep);
	    m_prescribedPressurePointerIndex = 0;
	}

	double* getPressurePointer();
	double* getPointerToFixedPressurePrescription();
	void setIsAtBoundary();
	bool isAtBoundary() const;
	double getPressure();
	void setPressure(const double pressure_in);
	void setPrescribedPressure(const double prescribedPressure);
	void setRestartPressureFromHistory();
	int getIndex() const;
	void setPythonControllerName(const std::string pythonParameterControllerName);
	bool hasUserDefinedExternalPythonScriptParameterController() const;
	std::string getPythonControllerName() const;

	double getHistoryPressure() const;
	void copyPressureToHistoryPressure();
	double getHistoryHistoryPressure() const;
	void copyHistoryPressureToHistoryHistoryPressure();
	bool hasHistoryPressure() const;
	void setHasHistoryPressure(const bool hasHistoryPressure);
	circuit_nodal_pressure_prescription_t getPressurePrescriptionType() const;
	void setPressurePrescriptionType(const circuit_nodal_pressure_prescription_t prescribedPressureType);
	int getPrescribedPressurePointerIndex() const;
	void setPrescribedPressurePointerIndex(const int prescribedPressurePointerIndex);
	double getFromPressureHistoryByTimestepIndex(const int timestepIndex) const;
	void appendToPressureHistory(const double pressure);



protected:
	double pressure;
	const int m_hstep;
private:
	double m_historyPressure;
	double m_historyHistoryPressure;
	bool m_hasHistoryPressure;
	circuit_nodal_pressure_prescription_t m_prescribedPressureType;
	int m_prescribedPressurePointerIndex;
	std::vector<double> m_entirePressureHistory;

	bool m_isAtBoundary;
	const int m_indexInInputData;
	bool m_hasPythonParameterController;
	std::string m_pythonParameterControllerName;
	double m_fixedPressure;
};

#endif