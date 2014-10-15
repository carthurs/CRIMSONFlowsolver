/*

 *
 *  Created on: Oct 7, 2014
 *      Author: klau, carthurs
 */

#include "multidom.h"
#include "fortranPointerManager.hpp"

// initialise the multidomain/LPN objects, this will need IFDEF for 3D and 1D codes
int boundaryCondition::bcCount = 0;
int RCR::numberOfInitialisedRCRs = 0;
int RCR::rcrtDatLinesReadSoFar = 0;

double boundaryCondition::getHop()
{
	return Hop;
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
      flowhist = new double [hstep];
      pressurehist = new double [hstep];
      
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
	flow_n = 0.0;
    flow_n1 = 0.0;
    pressure_n = 0.0;
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



void multidom_initialise(){
    
    

	std::string temp = "rcrt.dat";

	rcrtReader rcrtReader_instance(temp,1);
	rcrtReader_instance.readAndSplitMultiSurfaceInputFile();



	std::vector<std::pair<int,std::string>> surfaceList;
	// std::pair entry1;
	// entry1.first=3;
	// entry1.second="rcr";
    // std::pair <int,std::string> entry(3,"rcr"); 

	surfaceList.push_back(std::pair <int,std::string> (3,"rcr"));
	surfaceList.push_back(std::pair <int,std::string> (44,"netlist"));

	// Build a factory
	boundaryConditionFactory factory;

	std::vector<std::unique_ptr<boundaryCondition>> boundaryConditions;
	for (auto iterator=surfaceList.begin(); iterator !=surfaceList.end(); iterator++)
	{
		boundaryConditions.push_back(factory.createBoundaryCondition(iterator->first,iterator->second));
}

	// if (grcrbccom.numGRCRSrfs > 0)
	// {
	// 	for (int i=0; i < grcrbccom.numGRCRSrfs; i++)
	// 	{
	// 		std::cout << "writing here" << std::endl;
	// 		std::cout << grcrbccom.nsrflistGRCR[i+1] << std::endl;
	// 	}
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

