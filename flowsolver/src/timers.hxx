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
	const int m_initialTimestep;
	const int m_triggerTimestep;
	int m_currentTimestep;
};

#endif