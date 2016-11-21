#include "ImpedanceBoundaryCondition.hxx"
#include <iostream>
#include <boost/lexical_cast.hpp>
#include "fileReaders.hxx"

#define WHILE_LOOP_WARNING_ITERATION_COUNT 1000000

void checkSafetyCounter(const int safetyCounter, const std::stringstream& currentAction)
{
	if (safetyCounter > WHILE_LOOP_WARNING_ITERATION_COUNT){
		std::cout << "WARNING: Exceeded " << WHILE_LOOP_WARNING_ITERATION_COUNT << " line reads during: " << currentAction << "\n";
	}
}

std::pair<double,double> ImpedanceBoundaryCondition::computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt)
{
	//..try easiest convolution Q and Z constant per time step
	// do j=3,numTpoints+1
	//    ImpConvCoef(j,:) = ValueListImp(j-1,:)/numTpoints
	// enddo
	std::vector<double> ImpConvCoef(m_numberOfTimePointsInData + 1, 0.0);
	for (int index = 2; index < m_numberOfTimePointsInData; index++)
	{
		double interpolateToTime = (m_numberOfTimePointsInData + 1 - index) * delt; // Just mirroring the old Fortran here..
		ImpConvCoef.at(index) = m_impedanceLinearInterpolator.interpolateInTimeWithPeriodicExtrapolation(interpolateToTime) / (m_numberOfTimePointsInData - 1);
	}
      
	// ImpConvCoef(1,:) =zero
	ImpConvCoef.at(0) = 0.0;
	// ImpConvCoef(2,:) =zero
	ImpConvCoef.at(1) = 0.0;

	// ImpConvCoef(numTpoints+2,:) =  &
	//            ValueListImp(numTpoints+1,:)/numTpoints
	double interpolateToTime = 0.0;
	ImpConvCoef.at(m_numberOfTimePointsInData) = m_impedanceLinearInterpolator.interpolateInTimeWithPeriodicExtrapolation(interpolateToTime) / (m_numberOfTimePointsInData - 1);

// ! compensate for yalpha passed not y in Elmgmr()
      // ImpConvCoef(numTpoints+1,:)= ImpConvCoef(numTpoints+1,:) &
      //                   - ImpConvCoef(numTpoints+2,:)*(1.0-alfi)/alfi 
	ImpConvCoef.at(m_numberOfTimePointsInData - 1) -= ImpConvCoef.at(m_numberOfTimePointsInData) / m_generalizedAlphaMethodAlpha;
      // ImpConvCoef(numTpoints+2,:)= ImpConvCoef(numTpoints+2,:)/alfi
	ImpConvCoef.at(m_numberOfTimePointsInData) /= m_generalizedAlphaMethodAlpha;

       // formerly "ImpConvCoef(numTpoints+2,:)":
    // dp_dq = ValueListImp(numTpoints + 1) / numTpoints / m_generalizedAlphaMethodAlpha;
    double firstTimestepImpedance = m_timeVaryingImpedance.front().second;
    dp_dq = firstTimestepImpedance / (m_numberOfTimePointsInData - 1) / m_generalizedAlphaMethodAlpha; // very suspicious of this... but it is as the fortran code did it...

    // from the old subroutine pHist:
    double poldImp = 0.0;
    for (int ii=0; ii < m_numberOfTimePointsInData; ii++)
    {
        poldImp += m_qHistImp.at(ii) * ImpConvCoef.at(ii);
    }

    Hop = poldImp;
}

void ImpedanceBoundaryCondition::loadInputFiles()
{

	readImpedanceFlowHistory();
	readTimeDomainImpedance();
         
         // qHistImp.reserve(ntimeptpT+1);
         // qHistTry.reserve(ntimeptpT);
         // qHistTryF.reserve(ntimeptpT);
         // readTimeDomainImpedance(); //read impedance data and initialize begin/end values
         // do i=2,ntimeptpT
         //    call Impint((ntimeptpT-i+1)*Delt(1),i) //return Imp values in reverse order ZN->Z0 - fills ValueListImp - now m_timeVaryingImpedance
         // enddo

         // allocate (poldImp(0:MAXSURF)) //for pressure part that depends on the history only
}

void ImpedanceBoundaryCondition::readImpedanceFlowHistory()
{
	// from genini.f90 around line 100
	std::stringstream currentActionInfoMessage;
	currentActionInfoMessage << "Reading " << m_impedanceFlowHistoryFileName << std::endl;
	std::cout << currentActionInfoMessage;
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
	std::cout << currentActionInfoMessage;
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

	m_impedanceLinearInterpolator = LinearInterpolator(m_timeVaryingImpedance);
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