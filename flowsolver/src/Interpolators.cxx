#include "Interpolators.hxx"
#include <cmath>
#include <stdexcept>
#include <iostream>

// Linear interpolator, with constant value extrapolation for values outside the data range
// (returns the first value for time < first time point; returns the final value for time > final time point).
double LinearInterpolator::interpolateInTimeWithConstantExtrapolation(const double &timeToInterpolateTo) const
{
  // Linearly interpolates between pairs of (Time,Value) pairs, in time.
  // If we've reached a Time past the end of the last value in the array,
  // this just returns the final value of Value

  checkInterpolatorIsValid();

  if (m_dataToInterpolate.back().first <= timeToInterpolateTo)
  {
    // Case where we're off the final end of the time data series.
    return m_dataToInterpolate.back().second;
  }
  else if (m_dataToInterpolate.front().first >= timeToInterpolateTo)
  {
    // Case wehre we're at (or before) the beginning of the time data series.
    return m_dataToInterpolate.front().second;
  }

  return interpolate(timeToInterpolateTo);
}

// Linear interpolator, which repeats the data periodically if a timeToInterpolateTo outside of the
// data range is requested.
double LinearInterpolator::interpolateInTimeWithPeriodicExtrapolation(const double& timeToInterpolateTo) const
{
  checkInterpolatorIsValid();

  double dataStartTime = m_dataToInterpolate.front().first;
  double dataEndTime = m_dataToInterpolate.back().first;

  if (dataStartTime <= timeToInterpolateTo && timeToInterpolateTo < dataEndTime)
  {
    return interpolate(timeToInterpolateTo);
  }
  else if (timeToInterpolateTo >= dataEndTime) // timeToInterpolateTo > dataEndTime too if we got here
  {
    double dataWidthInTime = dataEndTime - dataStartTime;
    double timeToInterpolateTo_shiftedToDataInterval = std::fmod(timeToInterpolateTo - dataStartTime, dataWidthInTime) + dataStartTime;
    return interpolate(timeToInterpolateTo_shiftedToDataInterval);
  }
  else if (timeToInterpolateTo < dataStartTime)
  {
    double dataWidthInTime = dataEndTime - dataStartTime;
    double timeToInterpolateTo_shiftedToDataInterval = std::fmod(timeToInterpolateTo - dataStartTime, dataWidthInTime) + dataWidthInTime + dataStartTime;
    return interpolate(timeToInterpolateTo_shiftedToDataInterval);
  }
  else
  {
    throw std::logic_error("Logic failure in linear interpolator. Contact the developers.");
  }
}

double LinearInterpolator::interpolate(const double& timeToInterpolateTo) const
{
  checkInterpolatorIsValid();

  int timepointIndexAfterTimeToInterpolateTo = 0; // we're about to find this...
  while (m_dataToInterpolate.at(timepointIndexAfterTimeToInterpolateTo).first <= timeToInterpolateTo)
  {
    timepointIndexAfterTimeToInterpolateTo++;
  }

  double distanceThroughTimeInterval;
  distanceThroughTimeInterval = (timeToInterpolateTo - m_dataToInterpolate.at(timepointIndexAfterTimeToInterpolateTo-1).first) / (m_dataToInterpolate.at(timepointIndexAfterTimeToInterpolateTo).first - m_dataToInterpolate.at(timepointIndexAfterTimeToInterpolateTo-1).first);
  
  return m_dataToInterpolate.at(timepointIndexAfterTimeToInterpolateTo-1).second*(1.0 - distanceThroughTimeInterval) + m_dataToInterpolate.at(timepointIndexAfterTimeToInterpolateTo).second*distanceThroughTimeInterval;
}

void LinearInterpolator::checkInterpolatorIsValid() const
{
  if (m_dataToInterpolate.size() == 0)
  {
    throw std::out_of_range("LinearInterpolator not properly constructed before use: it has no data to interpolate.");
  }
}