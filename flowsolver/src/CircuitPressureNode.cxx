#include "CircuitPressureNode.hxx"

double* CircuitPressureNode::getPressurePointer()
{
	double* pressurePointer = &pressure;
	return pressurePointer;
}

double* CircuitPressureNode::getPointerToFixedPressurePrescription()
{
	// we should only be accessing this pointer for modification if it is a fixed-type pressure prescription
	assert(m_prescribedPressureType == Pressure_Fixed);
	double* pressurePointer = &m_fixedPressure;
	return pressurePointer;
}

bool CircuitPressureNode::hasUserDefinedExternalPythonScriptParameterController() const
{
	return m_hasPythonParameterController;
}

std::string CircuitPressureNode::getPythonControllerName() const
{
	assert(m_hasPythonParameterController);
	return m_pythonParameterControllerName;
}

// Sets the name "whateverName" of the nodal pressure controller to look for in the
// working directory: whateverName.py, containing class whateverName,
// with class method:
// newParamterValue = updateControl(self, oldParameterValue, delt).
//
// This should be a Python script.
void CircuitPressureNode::setPythonControllerName(const std::string pythonParameterControllerName)
{
	m_hasPythonParameterController = true;
	m_pythonParameterControllerName = pythonParameterControllerName;
}

int CircuitPressureNode::getIndex() const
{
	return m_indexInInputData;
}

void CircuitPressureNode::setIsAtBoundary()
{
	m_isAtBoundary = true;
}

bool CircuitPressureNode::isAtBoundary() const
{
	return m_isAtBoundary;
}

double CircuitPressureNode::getPressure()
{
	// If this is a prescribed fixed pressure, ensure we reset it to the original input value.
	// This has the additional benefit of stopping any drift in a supposedly-prescribed value.
	if (m_prescribedPressureType == Pressure_Fixed)
	{
		pressure = m_fixedPressure;
	}
	return pressure;
}

void CircuitPressureNode::setPressure(const double pressure_in)
{
	pressure = pressure_in;
}

void CircuitPressureNode::setPrescribedPressure(const double prescribedPressure)
{
	// We only do anything special with fixed-pressure values. There's nothing
	// to do here with other types of prescribed pressure, as they don't
	// remain fixed at a single value (so we needn't remember it in m_fixedPressure).
	if (m_prescribedPressureType == Pressure_Fixed)
	{
		m_fixedPressure = prescribedPressure;
	}
	else
	{
		pressure = prescribedPressure;
	}
}

void CircuitPressureNode::setRestartPressureFromHistory()
{
	pressure = m_entirePressureHistory.back();
}

void CircuitPressureNode::setHasHistoryPressure(const bool hasHistoryPressure)
{
	m_hasHistoryPressure = hasHistoryPressure;
}

void CircuitPressureNode::copyPressureToHistoryPressure()
{
	// std::cout << "setting history pressure to " << getPressure() << " from " << m_historyPressure << std::endl;
	m_historyPressure = getPressure();
}

void CircuitPressureNode::copyHistoryPressureToHistoryHistoryPressure()
{
	// used with the Kalman filter to store a history pressure from the step before the particle
	m_historyHistoryPressure = m_historyPressure;
}

double CircuitPressureNode::getHistoryPressure() const
{
	return m_historyPressure;
}

double CircuitPressureNode::getHistoryHistoryPressure() const
{
	// used with the Kalman filter to store a history pressure from the step before the particle
	return m_historyHistoryPressure;
}

bool CircuitPressureNode::hasHistoryPressure() const
{
	return m_hasHistoryPressure;
}

circuit_nodal_pressure_prescription_t CircuitPressureNode::getPressurePrescriptionType() const
{
	return m_prescribedPressureType;
}

void CircuitPressureNode::setPressurePrescriptionType(const circuit_nodal_pressure_prescription_t prescribedPressureType)
{
	m_prescribedPressureType = prescribedPressureType;
}

int CircuitPressureNode::getPrescribedPressurePointerIndex() const
{
	return m_prescribedPressurePointerIndex;
}

void CircuitPressureNode::setPrescribedPressurePointerIndex(const int prescribedPressurePointerIndex)
{
	m_prescribedPressurePointerIndex = prescribedPressurePointerIndex;
}

double CircuitPressureNode::getFromPressureHistoryByTimestepIndex(const int timestepIndex) const
{
	return m_entirePressureHistory.at(timestepIndex);
}

void CircuitPressureNode::appendToPressureHistory(const double pressure)
{
	// m_pressureHistoryBuffer.at(m_pressureHistoryBufferNextWriteIndex) = pressure;
	// m_pressureHistoryBufferNextWriteIndex++;
	// if (m_pressureHistoryBufferNextWriteIndex > m_pressureHistoryBufferSize - 1)
	// {
	// 	m_pressureHistoryBufferNextWriteIndex = 0;
	// 	m_entirePressureHistory.insert(m_entirePressureHistory.end(), m_pressureHistoryBuffer.begin(), m_pressureHistoryBuffer.end());
	// }

	m_entirePressureHistory.push_back(pressure);
}