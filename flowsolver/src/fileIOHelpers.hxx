#ifndef FILEIOHELPERS_HXX_
#define FILEIOHELPERS_HXX_
#include <vector>
#include <boost/shared_ptr.hpp>
#include "AbstractBoundaryCondition.hxx"
#include "ClosedLoopDownstreamSubsection.hxx"

void writeNetlistFlowsPressuresAndVolumes(const std::vector<boost::shared_ptr<AbstractBoundaryCondition>>& boundaryConditions, const std::vector<boost::shared_ptr<ClosedLoopDownstreamSubsection>> netlistDownstreamLoopClosingSubsections, int& nextTimestepWrite_start);
// void loadNetlistPressuresFlowsAndVolumesOnRestart(const std::vector<boost::shared_ptr<AbstractBoundaryCondition>>& boundaryConditions, const std::vector<boost::shared_ptr<ClosedLoopDownstreamSubsection>> netlistDownstreamLoopClosingSubsections, const int startingTimeStepIndex);

#endif