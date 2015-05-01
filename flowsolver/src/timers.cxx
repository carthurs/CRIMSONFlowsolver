#include "timers.hxx"
#include <iostream>

bool BasicTimer::hasTheTimeCome()
{
	bool theTimeHasCome = (m_currentTimestep >= m_triggerTimestep);
	std::cout << "Has the time come? : " << theTimeHasCome << std::endl;
	return (theTimeHasCome);
}

void BasicTimer::incrementTimer()
{
	std::cout << "timer says step is: " << m_currentTimestep << std::endl;
	m_currentTimestep++;
}