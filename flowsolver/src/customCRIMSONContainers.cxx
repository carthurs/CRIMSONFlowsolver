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

double ComponentParameterContainer::getInitialUnstressedVolume() const
{
	return m_initialUnstressedVolume;
}

void ComponentParameterContainer::setInitialUnstressedVolume(double inititalUnstressedVolume)
{
	m_initialUnstressedVolume = inititalUnstressedVolume;
}




void ComponentControlSpecificationContainer::addControlScript(const parameter_controller_t type, const std::string controlScriptName)
{
	m_parameterControllerTypesAndControlScriptNames.push_back(std::make_pair(type, controlScriptName));
}

int ComponentControlSpecificationContainer::getNumberOfControlScripts() const
{
	return m_parameterControllerTypesAndControlScriptNames.size();
}

parameter_controller_t ComponentControlSpecificationContainer::getControlTypeByIndexLocalToComponent(const int index) const
{
	return m_parameterControllerTypesAndControlScriptNames.at(index).first;
}

const std::string ComponentControlSpecificationContainer::getControlScriptNameByIndexLocalToComponent(const int index) const
{
	return m_parameterControllerTypesAndControlScriptNames.at(index).second;
}
