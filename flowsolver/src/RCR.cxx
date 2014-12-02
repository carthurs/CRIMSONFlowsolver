#include "RCR.hxx"
#include "fortranPointerManager.hxx"
#include "boundaryConditionManager.hxx"

// Statics
int RCR::numberOfInitialisedRCRs = 0;

double RCR::linInterpolateTimeData(const double &currentTime, const int timeDataLength)
{
  // Linearly interpolates between pairs of (Time,Value) pairs, in time.
  // If we've reached a Time past the end of the last value in the array,
  // this just returns the final value of Value

  if (timeDataPdist[timeDataLength-1].first <= currentTime)
  {
    // Case where we're off the final end of the time data series.
    return timeDataPdist[timeDataLength-1].second;
  }
  else if (timeDataPdist[0].first >= currentTime)
  {
    // Case wehre we're at (or before) the beginning of the time data series.
    return timeDataPdist[0].second;
  }

  int ii = 0;
  while (timeDataPdist[ii].first < currentTime)
  {
    ii++;
  }

  double distanceThroughTimeInterval;
  distanceThroughTimeInterval = (currentTime - timeDataPdist[ii-1].first) / (timeDataPdist[ii].first - timeDataPdist[ii-1].first);

  return timeDataPdist[ii-1].second*(1.0 - distanceThroughTimeInterval) + timeDataPdist[ii].second*distanceThroughTimeInterval;

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
    
    if (thisIsARestartedSimulation)
    {
      // Initialise the pressure using the value from the PHistRCR.dat.
      pressure_n = (boundaryConditionManager::Instance()->PHistReader)->getReadFileData(indexOfThisRCR+1,timdat.lstep);
    }
    else
    {
      pressure_n = *pressure_n_ptr;
    }

    // These are Hop and dp_dq now.
    // implicitcoeff = 0.0;
    // implicitcoeff_n1 = 0.0; 

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

// Here we step the actual discretised ODE for the RCR:
std::pair<double,double> RCR::computeImplicitCoefficients(int timestepNumber, double timeAtStepNplus1, double alfi_delt)
{

  double temp1;
  double temp2;
  std::pair<double,double> returnCoeffs;
  double timeAtStepN = delt*((double)timestepNumber);

  double rdn_1 = r2;
  double rp = r1;
  double compliance = c;

  double pdistn = linInterpolateTimeData(timeAtStepN,lengthOftimeDataPdist);
  double pdistn_1 = linInterpolateTimeData(timeAtStepNplus1,lengthOftimeDataPdist);

  // double pdistn = linInterpolateTimeData(timeDataPdist,timeAtStepN,lengthOftimeDataPdist);
  // double pdistn_1 = linInterpolateTimeData(timeDataPdist,timeAtStepNplus1,lengthOftimeDataPdist);

  // // parameters overwritten
  // // dirty hack for filtering
  //  rdn_1 = a%parameters_RCR(3,i)
  //  rp = a%parameters_RCR(1,i)
  //  compliance = a%parameters_RCR(2,i)

  //  pdistn = a%parameters_Pd
  //  pdistn_1 = a%parameters_Pd

  double denom = 1.0 + ((compliance*rdn_1)/alfi_delt);

  temp1 = rdn_1 + rp*(1.0 + ((compliance*rdn_1)/alfi_delt));

  temp2 = pressure_n + pdistn_1 - pdistn - rp*(*flow_n_ptr);
  temp2 = ((compliance*rdn_1)/alfi_delt)*temp2+ pdistn_1;

  returnCoeffs.first = temp1 / denom;
  returnCoeffs.second = temp2 / denom;
  
  return returnCoeffs;
}