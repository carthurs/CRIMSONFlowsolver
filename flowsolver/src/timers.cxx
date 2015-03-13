#include "timers.hxx"

bool BasicTimer::hasTheTimeCome()
{
	return (m_currentTimestep >= m_triggerTimestep);
}

void BasicTimer::incrementTimer()
{
	m_currentTimestep++;
}