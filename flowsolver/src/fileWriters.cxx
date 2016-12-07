#include <iomanip>
#include "fileWriters.hxx"
#include "mpi.h"

void BasicFileWriter::openFile()
{
	m_fileHandle = new std::ofstream(m_fileName.c_str(), std::ios::app);
}

void BasicFileWriter::writeStepIndex(int stepIndex)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank==0)
	{
		(*m_fileHandle) << std::setw(10) << std::left << stepIndex;
	}
}

void BasicFileWriter::writeToFile(double valueToWrite)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank==0)
	{
		(*m_fileHandle) << " " << std::setw(23) << std::left << valueToWrite;
	}
}

void BasicFileWriter::writeEndLine()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank==0)
	{
		(*m_fileHandle) << std::endl;
	}
}

void FileReplacingFileWriter::openFile()
{
	m_fileHandle = new std::ofstream(m_fileName.c_str(), std::ios::trunc);
}