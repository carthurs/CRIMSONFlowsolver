#ifndef RCR_HXX_
#define RCR_HXX_

#include "abstractBoundaryCondition.hxx"
#include "fileReaders.hxx"

class RCR : public abstractBoundaryCondition
{
public:
	RCR(int surfaceIndex_in)
	: abstractBoundaryCondition(surfaceIndex_in)
	{
		// Note the index of this RCR (zero-indexed), and count its existance
		// (the order of these two lines is correct!)
		indexOfThisRCR = numberOfInitialisedRCRs;
		numberOfInitialisedRCRs++;
		
		initialiseModel();
		rcrtReader* rcrtReader_instance = rcrtReader::Instance();
		r1 = rcrtReader_instance->getR1()[indexOfThisRCR];
		c = rcrtReader_instance->getC()[indexOfThisRCR];
		r2 = rcrtReader_instance->getR2()[indexOfThisRCR];
		timeDataPdist = rcrtReader_instance->getTimeDataPdist()[indexOfThisRCR];
		lengthOftimeDataPdist = rcrtReader_instance->getNumDataRCR()[indexOfThisRCR];
	}
	
 	std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timeAtStepNplus1, const double alfi_delt);

	
//  	procedure :: setimplicitcoeff_rcr => setimplicitcoeff_rcr
        // void setImplicitCoeff();
//     procedure :: updxvars_rcr => updxvars_rcr         
//     procedure :: writexvars_rcr => writexvars_rcr
//     procedure :: assign_ptrs_ext_rcr => assign_ptrs_ext_rcr

	~RCR()
	{
		numberOfInitialisedRCRs--;
	}
protected:
	double linInterpolateTimeData(const double &currentTime, const int timeDataLength);

private:
	void initialiseModel();
	static int numberOfInitialisedRCRs;
	int indexOfThisRCR;
	int pdmax;
	double r1; // Proximal resistance
	double c; // Capacitance (compliance)
	double r2; // Distal Resistance
	std::vector<std::pair<double,double>> timeDataPdist; // Time-varying disal pressure data
	int lengthOftimeDataPdist;
};

#endif