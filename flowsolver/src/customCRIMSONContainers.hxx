#ifndef CUSTOMCRIMSONNCONTAINERS_HXX_
#define CUSTOMCRIMSONNCONTAINERS_HXX_

class ComponentParameterContainer
{
public:
	ComponentParameterContainer()
	: m_hasInitialVolume(false)
	{
	}

	double getParameter() const;
	double getInitialVolume() const;

	void setParameter(double parameterValue);
	void setInitialVolume(double initialVolume);
private:
	double m_parameterValue;
	double m_initialVolume;
	bool m_hasInitialVolume;
};

#endif