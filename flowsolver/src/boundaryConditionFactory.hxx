#ifndef BOUNDARYCONDITIONFACTORY_HXX_
#define BOUNDARYCONDITIONFACTORY_HXX_

#include <boost/shared_ptr.hpp>
#include "abstractBoundaryCondition.hxx"

// Forward declarations:
class abstractBoundaryCondition;

class boundaryConditionFactory
 {
 public:
	static boost::shared_ptr<abstractBoundaryCondition> createBoundaryCondition(int surfaceIndex_in, std::string boundaryType);
 };

 #endif