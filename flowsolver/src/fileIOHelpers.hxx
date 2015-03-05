#ifndef FILEIOHELPERS_HXX_
#define FILEIOHELPERS_HXX_
#include <vector>
#include <boost/shared_ptr.hpp>
#include "abstractBoundaryCondition.hxx"

void writeNetlistFlowsPressuresAndVolumes(const std::vector<boost::shared_ptr<abstractBoundaryCondition>>& boundaryConditions, int& nextTimestepWrite_start);

#endif