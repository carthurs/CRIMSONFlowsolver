#ifndef TIMERS_HXX_
#define TIMERS_HXX_

class BasicTimer
{
public:
	BasicTimer(int initialTimestep, const int triggerTimestep)
	: m_initialTimestep(initialTimestep),
	m_triggerTimestep(triggerTimestep),
	m_currentTimestep(initialTimestep)
	{
	}
	
	bool hasTheTimeCome();
	void incrementTimer();
private:
	const int m_triggerTimestep;
	const int m_initialTimestep;
	int m_currentTimestep;
};

#endif