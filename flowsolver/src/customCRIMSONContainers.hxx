#ifndef CUSTOMCRIMSONNCONTAINERS_HXX_
#define CUSTOMCRIMSONNCONTAINERS_HXX_

#include <string>
#include <vector>
#include "datatypesInCpp.hxx"

class ComponentParameterContainer
{
public:
	ComponentParameterContainer()
	: m_hasInitialVolume(false)
	{
	}

	double getParameter() const;
	double getInitialVolume() const;
	double getInitialUnstressedVolume() const;

	void setParameter(double parameterValue);
	void setInitialVolume(double initialVolume);
	void setInitialUnstressedVolume(double inititalUnstressedVolume);
private:
	double m_parameterValue;
	double m_initialVolume;
	bool m_hasInitialVolume;
	double m_initialUnstressedVolume;
};

class ComponentControlSpecificationContainer
{
public:
	ComponentControlSpecificationContainer(){};
	void addControlScript(const parameter_controller_t type, const std::string controlScriptName);
	int  getNumberOfControlScripts() const;
	parameter_controller_t getControlTypeByIndexLocalToComponent(const int index) const;
	const std::string getControlScriptNameByIndexLocalToComponent(const int index) const;
private:
	std::vector<std::pair<parameter_controller_t, std::string>> m_parameterControllerTypesAndControlScriptNames;
};

#endif