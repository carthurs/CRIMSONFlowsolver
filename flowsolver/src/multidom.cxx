/*

 *
 *  Created on: Oct 7, 2014
 *      Author: klau, carthurs
 */

#include "common_c.h"
#include "multidom.h"
#include "fortranPointerManager.hxx"
#include <typeinfo>


// initialise the multidomain/LPN objects, this will need IFDEF for 3D and 1D codes
double boundaryCondition::getHop()
{
  return Hop;
}

double boundaryCondition::getdp_dq()
{
  return dp_dq;
}

 std::unique_ptr<boundaryCondition> boundaryConditionFactory::createBoundaryCondition (int surfaceIndex, std::string boundaryType)
 {

  if (boundaryType.compare("rcr") == int(0))
  {
    return std::unique_ptr<boundaryCondition> (new RCR(surfaceIndex));
  }
  else if (boundaryType.compare("netlist") == int(0))
    {
    return std::unique_ptr<boundaryCondition> (new netlist(surfaceIndex));
    }
  else
  {
    std::cout << "Unknown boundary type. Exiting.\n";
    std::exit(1);
  }

  // return returnObject;

}

double boundaryCondition::linInterpolateTimeData(const std::vector<std::pair<double,double>> &timeData, const double &currentTime, const int timeDataLength)
{
  // Linearly interpolates pairs of (Time,Value), in time.
  // If we've reached a Time past the end of the last value in the array,
  // this just returns the final value of Value

  if (timeData[timeDataLength].first <= currentTime)
  {
    // Case where we're off the final end of the time data series.
    return timeData[timeDataLength].second;
  }
  else if (timeData[0].first >= currentTime)
  {
    // Case wehre we're at (or before) the beginning of the time data series.
    return timeData[0].second;
  }

  int ii = 0;
  while (timeData[ii].first < currentTime)
  {
    ii++;
  }

  double distanceThroughTimeInterval;
  distanceThroughTimeInterval = (currentTime - timeData[ii-1].first) / (timeData[ii].first - timeData[ii-1].first);

  return timeData[ii-1].second*(1.0 - distanceThroughTimeInterval) + timeData[ii].second*distanceThroughTimeInterval;

}

void boundaryConditionManager::getImplicitCoeff_rcr(double* implicitCoeffs_toBeFilled)
{
  // This code is a bit tricky, becase FORTRAN/C++ interfacing doesn't yet support passing arrays which are sized
  // at run-time to C++ from FORTRAN. Therefore, I've had to just pass a pointer to the first entry, and then manage
  // dereferencing of that pointer manually to fill the whole array, but with the FORTRAN column-major array structure,
  // as opposed to the C++ row-major standard.
  int writeLocation = 0;
  for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(RCR))
    {
      implicitCoeffs_toBeFilled[writeLocation] = (*iterator)->getdp_dq();
      // +MAXSURF+1 here to move to the next column of the array (the +1 is annoying, and is because of weird design decisions in old FORTRAN code)
      implicitCoeffs_toBeFilled[writeLocation+MAXSURF+1] = (*iterator)->getHop();
      writeLocation++;
      // std::cout.precision(15);
      // std::cout << "set here in C++ hop: " << (*iterator)->getHop() << std::endl;
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCppGetImplicitCoeff_rcr(double*& implicitCoeffs_toBeFilled)
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->getImplicitCoeff_rcr(implicitCoeffs_toBeFilled);
}

void boundaryConditionManager::setSurfaceList(std::vector<std::pair<int,std::string>> surfaceList)
{
  // Build a factory
  boundaryConditionFactory factory;
  
  for (auto iterator=surfaceList.begin(); iterator !=surfaceList.end(); iterator++)
  {
    boundaryConditions.push_back(factory.createBoundaryCondition(iterator->first,iterator->second));
  }
    
  // std::cout << "the boundary condition has an r1 of: " << (*boundaryConditions[0]).tempDataTestFunction() << std::endl;
}

std::vector<std::unique_ptr<boundaryCondition>>* boundaryConditionManager::getBoundaryConditions()
{
    return &boundaryConditions;
}

void boundaryConditionManager::computeAllImplicitCoeff_solve(int timestepNumber)
{
  for (auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    (*iterator)->computeImplicitCoeff_solve(timestepNumber);
  }
}
// ---WRAPPED BY--->
extern "C" void callCppComputeAllImplicitCoeff_solve(int& timestepNumber)
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->computeAllImplicitCoeff_solve(timestepNumber);
}

void boundaryConditionManager::computeAllImplicitCoeff_update(int timestepNumber)
{
  for (auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    (*iterator)->computeImplicitCoeff_update(timestepNumber);
  }
}
// ---WRAPPED BY--->
extern "C" void callCppComputeAllImplicitCoeff_update(int& timestepNumber)
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->computeAllImplicitCoeff_update(timestepNumber);
}

//
/*  >>>>>>>>>>>>>>>>>>>>     RUUUULE BRITANNIAAAA! <<<<<<<<<<<<<<<<<<<<<<<<<<
ZZZ         888888888888888888888   ZZZZZZZZ  888888888888888888888   ZZZZZZ    
ZZZZZZ         888888888888888888   ZZZZZZZZ  888888888888888888   ZZZZZZ       
   ZZZZZZ         888888888888888   ZZZZZZZZ  888888888888888   ZZZZZZ         8
888   ZZZZZZ         888888888888   ZZZZZZZZ  888888888888   ZZZZZZ         8888
8888887  =ZZZZZ:        888888888   ZZZZZZZZ  888888888   ZZZZZZ        :8888888
8888888888   ZZZZZZ         88888   ZZZZZZZZ  88888   ZZZZZZ         88888888888
8888888888888   ZZZZZZ         88   ZZZZZZZZ  88   ZZZZZZ         88888888888888
8888888888888888   ZZZZZZ           ZZZZZZZZ    ZZZZZZ         88888888888888888
                                    ZZZZZZZZ                                    
                                    ZZZZZZZZ                                    
ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
                                    ZZZZZZZZ                                    
8888888888888888         ZZZZZZ     ZZZZZZZZ          ZZZZZZ   88888888888888888
8888888888888         ZZZZZZ   88   ZZZZZZZZ  88         ZZZZZZ   88888888888888
8888888888         ZZZZZZ   88888   ZZZZZZZZ  88888         ZZZZZZ   88888888888
8888888         ZZZZZ$  :88888888   ZZZZZZZZ  888888887        IZZZZZ.  88888888
888         ZZZZZZ   888888888888   ZZZZZZZZ  888888888888         ZZZZZZ   8888
         ZZZZZZ   888888888888888   ZZZZZZZZ  888888888888888         ZZZZZZ   8
      ZZZZZZ   888888888888888888   ZZZZZZZZ  888888888888888888         ZZZZZZ 
   ZZZZZZ   888888888888888888888   ZZZZZZZZ  888888888888888888888         ZZZZ
ZZZZZ$  :888888888888888888888888   ZZZZZZZZ  8888888888888888888888887        I
*/
void RCR::initialiseModel()
{


    
    
    // class(numericalrcr), intent(inout) :: this
    //   integer :: surfnum
    //   integer :: surflist(0:maxsurf)
    //   real*8 :: rcrcoeff(surfnum,3)
    //   integer :: pdmax, hstep
//       real*8 :: pdval(pdmax,2,surfnum) ! this array is initialised as zero
// !                                      ! then filled in, therefore values 
//                                       ! that are 0 then the index is > 2 
//                                       ! should be ignored  
    //   integer :: initrcr, ierr, fnum
    //   integer :: i, j ,k
    //   this%isactive = int(1)
    isactive = int(1);

    //   ! A label so we can identify this class (avoids a world of pain when working with a tower of derived types)
    //   this%classNameString = 'numericalrcr'  
// !     ! set surface numbers and lists   
//       this%surfnum = surfnum

    //   allocate(this%surfids(surfnum))
    //   this%surfids(1:surfnum) = surflist(1:surfnum)
    //   ! allocate history arrays
    //   if (lstep .gt. int(0)) then
    //      hstep = nstep + lstep
    //   else
    //      hstep = nstep
    //   end if      
    // allocate(this%flowhist(hstep+1,surfnum))       
    // allocate(this%pressurehist(hstep+1,surfnum))             


    //   ! zero history arrays
    //   this%flowhist(:,:) = real(0.0,8)
      // flowhist = new double [hstep];
      // pressurehist = new double [hstep];
      
      for(int ii=0; ii<hstep; ii++)
      {
          flowhist[ii] = 0.0;
          pressurehist[ii] = 0.0;
      }

    //   this%pressurehist(:,:) = real(0.0,8)
    //   ! set flow and pressure file names

    //   write(this%flowfile,'(a)') 'QHistRCR.dat'
    flowfile = "QHistRCR.dat";

    //   write(this%pressurefile,'(a)') 'PHistRCR.dat'      
    pressurefile = "PHistRCR.dat";

    //   ! allocate arrays for input parameters & data
    //   allocate(this%rcrparams(surfnum))
    //   allocate(this%parameters_RCR(3,surfnum)) ! testing
    //   ! allocate and zero other arrays 
    // allocate(this%surfarea(surfnum))
    // allocate(this%flow_n(surfnum))
    // allocate(this%flow_n1(surfnum)) !! flow at n+1 / n+alf NOT USED
    // allocate(this%pressure_n(surfnum))      
    // allocate(this%implicitcoeff(surfnum,2)) 
    // allocate(this%implicitcoeff_n1(surfnum,2)) 

    surfarea = 0.0;

    // here we set the initial values of the flow and pressure using the pointers to the multidomaincontainer.
    // NB: Need to add a method in fortran to set a value for non-zero restarting!
    flow_n_ptr = fortranBoundaryDataPointerManager::Get()->boundaryFlows.at(surfaceIndex);
    pressure_n_ptr = fortranBoundaryDataPointerManager::Get()->boundaryPressures.at(surfaceIndex);
    
    std::cout << "just set pressure and flow: " << *pressure_n_ptr << " " << *flow_n_ptr << std::endl;

    flow_n1 = 0.0;
    
    
    implicitcoeff = 0.0;
    implicitcoeff_n1 = 0.0; 

    // ! initialise reservoir pressure  
    // if (initrcr) then
      //    ! set int
      //    this%init_pRes = 1
      //    allocate(this%pRes_n(surfnum))
      //    ! open rcr.x.dat
      //    open(fnum, file='rcrt.x.dat', status='old', iostat=ierr)         
      //    do i = 1, surfnum
      //       read(fnum,*) this%pRes_n(i)
      //    end do
      //    close(fnum)
      // else 
      //    this%init_pRes = 0
      // end if

      // ! zero variables 
      // this%surfarea(:) = real(0.0,8)
      // this%flow_n(:) = real(0.0,8)
      // this%flow_n1(:) = real(0.0,8)
      // this%pressure_n(:) = real(0.0,8)
      // this%implicitcoeff(:,:) = real(0.0,8)
      // this%implicitcoeff_n1(:,:) = real(0.0,8)



      // do i = 1, surfnum
      //    ! set rcr parameters
      //    this%rcrparams(i)%rp = rcrcoeff(i,1)
      //    this%rcrparams(i)%c = rcrcoeff(i,2)
      //    this%rcrparams(i)%rd = rcrcoeff(i,3)
      //    this%parameters_RCR(3,i) = rcrcoeff(i,3)
      //    this%parameters_RCR(1,i) = rcrcoeff(i,1)
      //    this%parameters_RCR(2,i) = rcrcoeff(i,2)
      //    ! count time points 
      //    k = int(1)
      //    do j = 2, pdmax
      //      if (pdval(j,1,i) .lt. pdval(j-1,1,i)) then
      //         exit
      //      else
      //         k = k + 1
      //      end if
      //    end do        
      //    ! allocate i'th entry with k x 2 size
      //    allocate(this%rcrparams(i)%pd%v(k,2))     
      //    ! set values
      //    do j = 1, k
      //       this%rcrparams(i)%pd%v(j,1) = pdval(j,1,i)
      //       this%rcrparams(i)%pd%v(j,2) = pdval(j,2,i)
      //    end do
      // end do

}

void RCR::computeImplicitCoeff_solve(int timestepNumber)
{
  std::pair<double,double> temp;

  double timeAtStepNplus1 = delt*((double)timestepNumber+timdat.alfi);
  double alfi_delt = timdat.alfi*delt;

  temp = computeImplicitCoefficients(timestepNumber, timeAtStepNplus1, alfi_delt);

  dp_dq = temp.first;
  Hop = temp.second;
  std::cout << "dp_dq in C++: " << dp_dq << std::endl;
  std::cout << "Hop in C++: " << Hop << std::endl;
}

void RCR::computeImplicitCoeff_update(int timestepNumber)
{
  std::pair<double,double> temp;

  double timeAtStepNplus1 = delt*((double)timestepNumber+1.0);
  double alfi_delt = delt;

  temp = computeImplicitCoefficients(timestepNumber, timeAtStepNplus1, alfi_delt);

  dp_dq_n1 = temp.first;
  Hop_n1 = temp.second;

  std::cout << "dp_dq_n1 in C++: " << dp_dq_n1 << std::endl;
  std::cout << "Hop_n1 in C++: " << Hop_n1 << std::endl;
}

std::pair<double,double> RCR::computeImplicitCoefficients(int timestepNumber, double timeAtStepNplus1, double alfi_delt)
{
  double temp1;
  double temp2;
  std::pair<double,double> returnCoeffs;
  double timeAtStepN = delt*((double)timestepNumber);

  double rdn_1 = r2;
  double rp = r1;
  double compliance = c;

  double pdistn = linInterpolateTimeData(timeDataPdist,timeAtStepN,lengthOftimeDataPdist);
  double pdistn_1 = linInterpolateTimeData(timeDataPdist,timeAtStepNplus1,lengthOftimeDataPdist);

  // // parameters overwritten
  // // dirty hack for filtering
  //  rdn_1 = a%parameters_RCR(3,i)
  //  rp = a%parameters_RCR(1,i)
  //  compliance = a%parameters_RCR(2,i)

  //  pdistn = a%parameters_Pd
  //  pdistn_1 = a%parameters_Pd

  double denom = 1.0 + ((compliance*rdn_1)/alfi_delt);

  temp1 = rdn_1 + rp*(1.0 + ((compliance*rdn_1)/alfi_delt));

  temp2 = (*pressure_n_ptr) + pdistn_1 - pdistn - rp*(*flow_n_ptr);
  temp2 = ((compliance*rdn_1)/alfi_delt)*temp2+ pdistn_1;

  returnCoeffs.first = temp1 / denom;
  returnCoeffs.second = temp2 / denom;
  
  return returnCoeffs;
}

int boundaryCondition::bcCount = 0;
int RCR::numberOfInitialisedRCRs = 0;

void multidom_initialise(){

  rcrtReader* rcrtReader_instance = rcrtReader::Instance();
  rcrtReader_instance->setFileName("rcrt.dat");
  rcrtReader_instance->readAndSplitMultiSurfaceInputFile();

  std::vector<std::pair<int,std::string>> surfaceList;

  // loop through rcr boundaries listed in the input file, surface numbers read from the common_c.h
  for (int i = 0; i < grcrbccom.numGRCRSrfs; i++)
  {
    surfaceList.push_back(std::pair <int,std::string> (grcrbccom.nsrflistGRCR[i+1],"rcr"));
  }
  // Write loops here for all the other surface types!
  
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->setSurfaceList(surfaceList);
  
  std::vector<std::unique_ptr<boundaryCondition>>* retrievedBoundaryConditions;
  retrievedBoundaryConditions = boundaryConditionManager_instance->getBoundaryConditions();
  

  // surfaceList.push_back(std::pair <int,std::string> (3,"rcr"));
  // surfaceList.push_back(std::pair <int,std::string> (4,"rcr"));
  // surfaceList.push_back(std::pair <int,std::string> (44,"netlist")); 


   std::cout << "the boundary condition has an r1 of: " << (*retrievedBoundaryConditions)[0]->tempDataTestFunction() << std::endl;

  // if (grcrbccom.numGRCRSrfs > 0)
  // {
  //  for (int i=0; i < grcrbccom.numGRCRSrfs; i++)
  //  {
  //    std::cout << "writing here" << std::endl;
  //    std::cout << grcrbccom.nsrflistGRCR[i+1] << std::endl;
  //  }
  // }

}


// set pointer to fortran arrays
void multidom_iter_initialise(){

}

void multidom_iter_step(){

}

void multidom_iter_finalise(){

}

void multidom_finalise(){

}

