#ifndef _PARTITION_H_
#define _PARTITION_H_

#include <vector>

extern void Partition_Problem( int );

void writeNodalOutputIndicesForProcessor(const int processorId, const std::vector<std::vector<int>>& processorZeroIndexToLocalNodalFlowAndPressureOutputNodeIndices_in);

#endif
