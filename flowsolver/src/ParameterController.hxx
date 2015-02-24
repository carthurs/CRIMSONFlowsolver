#ifndef PARAMETERCONTROLLER_HXX_
#define PARAMETERCONTROLLER_HXX_

#include <cmath>

// This class can be connected to a netlist LPN parameter
// (e.g. resistance for a component, or compliance for a compliance chamber).
// It will adjust that parameter according to necessary control inputs.
class AbstractParameterController
{
public:
	AbstractParameterController(double* const parameterToControl)
	: mp_parameterToControl(parameterToControl), // Set the pointer to the parameter we will control
	m_originalValueOfParameter(*parameterToControl) // Save the original state of the parameter
	{
	}
	virtual ~AbstractParameterController()
	{
	}

	// Function to adjust the controlled parameter as necessary.
	virtual void updateControl() = 0;
protected:
	double* const mp_parameterToControl;
	const double m_originalValueOfParameter;
};


class LeftVentricularElastanceController : public AbstractParameterController
{
public:
	LeftVentricularElastanceController(double* const parameterToControl, const double delt)
	: AbstractParameterController(parameterToControl),
	m_delt(delt)
	{
		m_periodicTime = 0.0; //\todo think about this for restarts!
		m_timeToMaximumElastance = 0.2782;
		m_timeToRelax = 0.1391;
		m_minimumElastance = 41.0246;
		m_maximumElastance = 3085.6;
		m_heartPeriod = 0.86;
	}

	void updateControl()
	{
		updatePeriodicTime();
		// adjust the controlled elastance:
		*mp_parameterToControl = 1.0/getElastance();
	}
private:
	const double m_delt;
	
	double m_periodicTime;
	double m_heartPeriod;
	
	double m_timeToMaximumElastance;
	double m_timeToRelax;

	double m_minimumElastance;
	double m_maximumElastance;

	double getElastance();

	void updatePeriodicTime();
};

#endif