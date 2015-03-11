#ifndef FILEWRITERS_HXX_
#define FILEWRITERS_HXX_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include "multidom.hxx"


class basicFileWriter
{
public:
	basicFileWriter()
	{
	}
	
	void setFileName(std::string fileNameIn)
	{
	    fileName = fileNameIn;
		fileHandle = new std::ofstream(fileName.c_str(), std::ios::app);

		fileHandle->precision(15);
		
		if (!fileHandle->is_open())
		{
			std::cout << "Failed to open " << fileName << " for writing!" << std::endl;
			std::exit(1);
		}
	}

	void writeStepIndex(int stepIndex);
	void writeToFile(double valueToWrite);
	void writeEndLine();
	
	~basicFileWriter()
	{
		fileHandle->close();
		delete fileHandle;
	}
protected:
	std::ofstream* fileHandle;
	std::string currentLine;
	std::string fileName;
	// double* arrayToWriteOut;

	int numColumns;
	bool appendLine();
private:
};

#endif
