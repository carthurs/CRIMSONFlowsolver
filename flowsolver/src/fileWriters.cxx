#include <iomanip>
#include "fileWriters.hxx"
#include "mpi.h"


void basicFileWriter::writeStepIndex(int stepIndex)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank==0)
	{
		(*fileHandle) << std::setw(10) << std::left << stepIndex;
	}
}

void basicFileWriter::writeToFile(double valueToWrite)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank==0)
	{
		(*fileHandle) << " " << std::setw(23) << std::left << valueToWrite;
	}
}

void basicFileWriter::writeEndLine()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank==0)
	{
		(*fileHandle) << std::endl;
	}
}