#ifndef FILEWRITERS_HXX_
#define FILEWRITERS_HXX_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include "multidom.hxx"


class BasicFileWriter
{
public:
	BasicFileWriter()
	{
	}
	
	void setFileName(std::string fileNameIn)
	{
	    m_fileName = fileNameIn;
		openFile();

		m_fileHandle->precision(15);
		
		if (!m_fileHandle->is_open())
		{
			std::cout << "Failed to open " << m_fileName << " for writing!" << std::endl;
			std::exit(1);
		}
	}

	void writeStepIndex(int stepIndex);
	void writeToFile(double valueToWrite);
	void writeEndLine();
	
	virtual ~BasicFileWriter()
	{
		m_fileHandle->close();
		delete m_fileHandle;
	}
protected:
	std::ofstream* m_fileHandle;
	std::string m_fileName;
	// double* arrayToWriteOut;

	int numColumns;
	bool appendLine();
private:
	virtual void openFile();
};

class FileReplacingFileWriter : public BasicFileWriter
{
private:
	void openFile() override;
};

#endif
