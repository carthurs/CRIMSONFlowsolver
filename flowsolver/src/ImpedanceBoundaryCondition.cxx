#include "ImpedanceBoundaryCondition.hxx"
#include <iostream>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include "fileReaders.hxx"
#include "fileWriters.hxx"

#define WHILE_LOOP_WARNING_ITERATION_COUNT 1000000

void checkSafetyCounter(const int safetyCounter, const std::stringstream& currentAction)
{
	if (safetyCounter > WHILE_LOOP_WARNING_ITERATION_COUNT){
		std::cout << "WARNING: Exceeded " << WHILE_LOOP_WARNING_ITERATION_COUNT << " line reads during: " << currentAction.str() << "\n";
	}
}

std::pair<double,double> ImpedanceBoundaryCondition::computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt)
{
    double firstTimestepImpedance = m_timeVaryingImpedance.front().second;

    // This is just the first entry of the impedance; it is really part
    // of the convolution, but the flow value it needs to be multiplied
    // by in the convolution is not known yet; the 3D domain will
    // find that value and multiply it by this dp_dq during the solver.
    //
    // No delt required in this as it's a discrete convolution
    dp_dq = firstTimestepImpedance / m_generalizedAlphaMethodAlpha;

    // c.f. the old subroutine pHist:
    double poldImp = 0.0;
    for (size_t ii = 0; ii < m_numberOfTimePointsInData - 1; ii++)
    {
        double convolutionTime = (ii+1) * delt;
        // No delt required in this as it's a discrete convolution
        poldImp += m_qHistImp.at(m_numberOfTimePointsInData - ii - 1) * mp_impedanceLinearInterpolator->interpolateInTimeWithPeriodicExtrapolation(convolutionTime) / m_generalizedAlphaMethodAlpha;
    }

    Hop = poldImp;
    return std::make_pair(dp_dq, Hop);
}

void ImpedanceBoundaryCondition::loadInputFiles()
{

	readImpedanceFlowHistory();
	readTimeDomainImpedance();
}

void ImpedanceBoundaryCondition::readImpedanceFlowHistory()
{
	// from genini.f90 around line 100
	std::stringstream currentActionInfoMessage;
	currentActionInfoMessage << "Reading " << m_impedanceFlowHistoryFileName << std::endl;
	std::cout << currentActionInfoMessage.str();

	SimpleFileReader simpleFileReader(m_impedanceFlowHistoryFileName);

	// Read the file into m_qHistImp:
	int safetyCounter = 0;
	while (true)
	{
		bool fileNotEnded = true;
		std::string valueJustRead = simpleFileReader.getNextDataSplitBySpacesOrEndOfLine(fileNotEnded);
		if (!fileNotEnded){
			break;
		}
		m_qHistImp.push_back( boost::lexical_cast<double>(valueJustRead) );

		safetyCounter++;
		checkSafetyCounter(safetyCounter, currentActionInfoMessage);
	}

	m_numberOfTimePointsInData = m_qHistImp.size(); // note that this is actually equal to ntimeptpT+1 in genini.
}

// -----------------------------------------------------------------------
//    initialize the impedance boundary condition:
//    read the data in initImpt
//    interpolate the data to match the process time step in Impint
// -----------------------------------------------------------------------

// subroutine initImpt()
void ImpedanceBoundaryCondition::readTimeDomainImpedance()
{
	std::stringstream currentActionInfoMessage;
	currentActionInfoMessage << "Reading " << m_timeDomainImpedanceFileName << std::endl;
	std::cout << currentActionInfoMessage.str();

	SimpleFileReader timeDomainImpedanceDatReader(m_timeDomainImpedanceFileName);

	int safetyCounter = 0;
	while (true)
	{
		bool fileNotEnded = true;

		// Read the time value from the first column of the file
		std::string timeValueJustReadFromFileColumnOne = timeDomainImpedanceDatReader.getNextDataSplitBySpacesOrEndOfLine(fileNotEnded);
		if (!fileNotEnded){
			break;
		}
		
		// Read the corresponding impedance value from teh second column of the file:
		std::string impedanceValueJustReadFromFileColumnTwo = timeDomainImpedanceDatReader.getNextDataSplitBySpacesOrEndOfLine(fileNotEnded);
		if (!fileNotEnded)
		{
			// A failure here is a disaster, so throw. There was a time value without a corresponding impedance value at the end of the file.
			std::stringstream errorMessage;
			errorMessage << "File " << m_timeDomainImpedanceFileName << "ended unexpectedly." << std::endl;
			throw std::runtime_error(errorMessage.str());
		}
		double timeValue = boost::lexical_cast<double>(timeValueJustReadFromFileColumnOne);
		double impedanceValue = boost::lexical_cast<double>(impedanceValueJustReadFromFileColumnTwo);

		// Fill the data structure formerly known as ValueImpt:
		m_timeVaryingImpedance.push_back( std::make_pair(timeValue, impedanceValue) );

		safetyCounter++;
		checkSafetyCounter(safetyCounter, currentActionInfoMessage);
	}

	#ifdef BOOST_NO_CXX11_SMART_PTR
	mp_impedanceLinearInterpolator = boost::shared_ptr<LinearInterpolator>(new LinearInterpolator(m_timeVaryingImpedance));
	#else
	// prefer unique_ptr if it's available
	mp_impedanceLinearInterpolator = std::unique_ptr<LinearInterpolator>(new LinearInterpolator(m_timeVaryingImpedance));
	#endif
}

// Replaces the code in subroutine UpdHistConv
void ImpedanceBoundaryCondition::updateStoredFlowHistory()
{
	m_qHistImp.pop_front();
	m_qHistImp.front() = 0.0;
	m_qHistImp.push_back(*flow_n_ptrs.at(0));
}

void ImpedanceBoundaryCondition::writeImpedanceSpecificFlowHistory() const
{
	FileReplacingFileWriter fileWriter;
	fileWriter.setFileName(m_impedanceFlowHistoryFileName);
	for (double value : m_qHistImp)
	{
		fileWriter.writeToFile(value);
		fileWriter.writeEndLine();
	}
}

void ImpedanceBoundaryCondition::finaliseAtEndOfTimestep()
{
	updateStoredFlowHistory();
	writeImpedanceSpecificFlowHistory();
}
