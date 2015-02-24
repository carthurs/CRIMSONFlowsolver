#include <iomanip>
#include "fileWriters.hxx"


void basicFileWriter::writeStepIndex(int stepIndex)
{
	(*fileHandle) << std::setw(10) << std::left << stepIndex;
}

void basicFileWriter::writeToFile(double valueToWrite)
{
	(*fileHandle) << " " << std::setw(23) << std::left << valueToWrite;
}

void basicFileWriter::writeEndLine()
{
	(*fileHandle) << std::endl;
}