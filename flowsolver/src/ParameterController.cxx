#include "ParameterController.hxx"
#include <iostream>

void LeftVentricularElastanceController::updateControl()
{
	updatePeriodicTime();
	// adjust the controlled elastance:
	*mp_parameterToControl = getElastance();
}

double LeftVentricularElastanceController::getElastance()
{
	// *** analytical elastance function from:
	//     pope, s. r.; ellwein, l. m.; zapata, c. l.; novak, v.; kelley, c. t. & olufsen, m. s.  
	//     estimation and identification of parameters in a lumped cerebrovascular model.
	//     math biosci eng, 2009, 6, 93-115
	
	double elastance;

	// This is the elastance function. It's defined piecewise:
	if ( m_periodicTime <= m_timeToMaximumElastance )
	{
 		elastance = m_minimumElastance
            + 0.5*(m_maximumElastance - m_minimumElastance)
            * (1.0 - cos((m_periodicTime*M_PI)/m_timeToMaximumElastance));
    }
  	else if ( m_periodicTime <= (m_timeToMaximumElastance + m_timeToRelax) )
  	{
     	elastance = m_minimumElastance
            + 0.5*(m_maximumElastance-m_minimumElastance)
            * (1.0 + cos((m_periodicTime-m_timeToMaximumElastance)*(M_PI/m_timeToRelax)));
    }
	else if ( m_periodicTime > (m_timeToMaximumElastance + m_timeToRelax) )
	{
		elastance = m_minimumElastance;
	}

	std::cout << "elastance was: "<< elastance << std::endl;

  return elastance;
}

void LeftVentricularElastanceController::updatePeriodicTime()
{
	m_periodicTime = m_periodicTime + m_delt;
	// Keep m_periodicTime in the range [0,m_heartPeriod) :
	if (m_periodicTime >= m_heartPeriod)
	{
		m_periodicTime = m_periodicTime - m_heartPeriod;
	}
	std::cout << "m_periodicTime was: "<< m_periodicTime << std::endl;
}

void BleedController::updateControl()
{
	bool m_bleedingOn = mp_timer->hasTheTimeCome();
	if (m_bleedingOn)
	{
		*mp_parameterToControl = 0.001; // set the resistance / compliance to be tiny (depending on the type of component we're controlling here...)
	}
	mp_timer->incrementTimer();
}