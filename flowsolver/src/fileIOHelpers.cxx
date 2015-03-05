#include "fileIOHelpers.hxx"
#include "NetlistBoundaryCondition.hxx"
#include "fileWriters.hxx"
#include <sstream>

void writeNetlistFlowsPressuresAndVolumes(const std::vector<boost::shared_ptr<abstractBoundaryCondition>>& boundaryConditions, int& nextTimestepWrite_start)
{
  int nextTimestepWrite_end;
  for (auto boundaryCondition=boundaryConditions.begin(); boundaryCondition!=boundaryConditions.end(); boundaryCondition++)
  {
    NetlistBoundaryCondition* netlistBoundaryCondition = dynamic_cast<NetlistBoundaryCondition*>(boundaryCondition->get());
    if (netlistBoundaryCondition != NULL) // if the boundaryCondition was a NetlistBoundaryCondition (or subclass thereof)
    {
      // All the following writes can use the same nextTimestepWrite_end to determine how far to go when looking for data to write to the file:
      nextTimestepWrite_end = netlistBoundaryCondition->getCircuitDescription()->components.at(0)->m_entireFlowHistory.size();
      
      {
        // Write the netlistPressures_surface_X.dat
        basicFileWriter boundaryConditionPressureHistoryWriter;
        std::stringstream filenameForThisBoundary;
        filenameForThisBoundary << "netlistPressures_surface_" << (*boundaryCondition)->getSurfaceIndex() << ".dat";
        boundaryConditionPressureHistoryWriter.setFileName(filenameForThisBoundary.str());

        for (int stepToWrite=nextTimestepWrite_start; stepToWrite<nextTimestepWrite_end; stepToWrite++)
        {
          boundaryConditionPressureHistoryWriter.writeStepIndex(stepToWrite);
          for (auto node=netlistBoundaryCondition->getCircuitDescription()->mapOfPressureNodes.begin(); node!=netlistBoundaryCondition->getCircuitDescription()->mapOfPressureNodes.end(); node++)
          {
            boundaryConditionPressureHistoryWriter.writeToFile(node->second->m_entirePressureHistory.at(stepToWrite));
          }
          boundaryConditionPressureHistoryWriter.writeEndLine();
        }
      }

      {
        // Write the netlistFlows_surface_X.dat
        basicFileWriter boundaryConditionFlowHistoryWriter;
        std::stringstream filenameForThisBoundary;
        filenameForThisBoundary << "netlistFlows_surface_" << (*boundaryCondition)->getSurfaceIndex() << ".dat";
        boundaryConditionFlowHistoryWriter.setFileName(filenameForThisBoundary.str());

        for (int stepToWrite=nextTimestepWrite_start; stepToWrite<nextTimestepWrite_end; stepToWrite++)
        {
          boundaryConditionFlowHistoryWriter.writeStepIndex(stepToWrite);
          for (auto component=netlistBoundaryCondition->getCircuitDescription()->components.begin(); component!=netlistBoundaryCondition->getCircuitDescription()->components.end(); component++)
          {
            boundaryConditionFlowHistoryWriter.writeToFile((*component)->m_entireFlowHistory.at(stepToWrite));
          }
          boundaryConditionFlowHistoryWriter.writeEndLine();
        }
      }


      {
        // Write the volumes of the pressure chambers, as netlistVolumes_surface_X.dat
        basicFileWriter boundaryConditionVolumeHistoryWriter;
        std::stringstream filenameForThisBoundary;
        filenameForThisBoundary << "netlistVolumes_surface_" << (*boundaryCondition)->getSurfaceIndex() << ".dat";
        boundaryConditionVolumeHistoryWriter.setFileName(filenameForThisBoundary.str());

        for (int stepToWrite=nextTimestepWrite_start; stepToWrite<nextTimestepWrite_end; stepToWrite++)
        {
          boundaryConditionVolumeHistoryWriter.writeStepIndex(stepToWrite);
          for (auto component=netlistBoundaryCondition->getCircuitDescription()->mapOfComponents.begin(); component!=netlistBoundaryCondition->getCircuitDescription()->mapOfComponents.end(); component++)
          {
            VolumeTrackingPressureChamber* pressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (component->second.get());
            // If this component is actually a volume chamber, so it actually has a volume history we can write to the file:
            if (pressureChamber != NULL)
            {
              boundaryConditionVolumeHistoryWriter.writeToFile(pressureChamber->m_entireVolumeHistory.at(stepToWrite));
            }
          }
          boundaryConditionVolumeHistoryWriter.writeEndLine();
        }
      }


    }
  }
  nextTimestepWrite_start = nextTimestepWrite_end;
}