#ifndef FILEWRITERS_HXX_
#define FILEWRITERS_HXX_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include "multidom.h"


class basicFileWriter
{
public:
	basicFileWriter()
	{
	}
	
	void setFileName(std::string fileNameIn)
	{
	    fileName = fileNameIn;
		fileHandle = new std::ofstream(fileName, std::ios::app);
		
		if (!fileHandle->is_open())
		{
			std::cout << "Failed to open for writing " << fileName << "!" << std::endl;
			std::exit(1);
		}
	}

	// void setArrayToWriteOut(double* arrayName)
	// {
	// 	arrayToWriteOut = arrayName;
	// }

	void writeHistories_rcr(std::vector<std::unique_ptr<boundaryCondition>> &boundaryConditions);

	// void setNumColumns(int numberOfColumns)
	// {
	// 	numColumns = numberOfColumns;
	// }
	
	~basicFileWriter()
	{
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