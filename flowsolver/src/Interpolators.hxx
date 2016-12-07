#ifndef INTERPOLATORS_HXX_
#define INTERPOLATORS_HXX_

#include "datatypesInCpp.hxx"

class LinearInterpolator
{
public:
	LinearInterpolator(const TimeValuePairVector dataToInterpolate)
	: m_dataToInterpolate(dataToInterpolate)
	{
	}

	LinearInterpolator()
	{
	}
	
	double interpolateInTimeWithConstantExtrapolation(const double &timeToInterpolateTo) const;
	double interpolateInTimeWithPeriodicExtrapolation(const double &timeToInterpolateTo) const;
private:
	const TimeValuePairVector m_dataToInterpolate;

	double interpolate(const double& timeToInterpolateTo) const;
	inline void checkInterpolatorIsValid() const;
};

#endif