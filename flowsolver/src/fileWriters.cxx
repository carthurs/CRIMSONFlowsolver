#include <iomanip>
#include "fileWriters.hxx"
#include "mpi.h"


void BasicFileWriter::writeStepIndex(int stepIndex)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank==0)
	{
		(*fileHandle) << std::setw(10) << std::left << stepIndex;
	}
}

void BasicFileWriter::writeToFile(double valueToWrite)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank==0)
	{
		(*fileHandle) << " " << std::setw(23) << std::left << valueToWrite;
	}
}

void BasicFileWriter::writeEndLine()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank==0)
	{
		(*fileHandle) << std::endl;
	}
}