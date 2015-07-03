#include "customCRIMSONContainers.hxx"
#include <stdexcept>

double ComponentParameterContainer::getParameter() const
{
	return m_parameterValue;
}

double ComponentParameterContainer::getInitialVolume() const
{
	if (!m_hasInitialVolume)
	{
		throw std::runtime_error("EE: Attempted to access an initial stored volume for a component which either doesn't store volume, or has not been given an initial volume by the inpt data.");
	}
	return m_initialVolume;
}

void ComponentParameterContainer::setParameter(double parameterValue)
{
	m_parameterValue = parameterValue;
}

void ComponentParameterContainer::setInitialVolume(double initialVolume)
{
	m_hasInitialVolume = true;
	m_initialVolume = initialVolume;
}