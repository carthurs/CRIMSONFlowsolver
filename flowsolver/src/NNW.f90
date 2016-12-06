module nnw
! module containing all the subroutines related to; 
! non newtonian flow SL, MAF December '16

implicit none

private ! everything is private by default


public :: initialize_NNW

integer, public :: typeViscConstModel
integer, public :: nparamViscConstModel
real*8, allocatable, public :: paramViscConstModel(:)

contains

subroutine initialize_NNW(nnwType)
implicit none
integer :: nnwType
integer :: ierr, rerr
integer :: fnum = 145
integer :: nparam, i
character(len=100) :: filename
logical :: file_exists

if (nnwType.eq.3) then
    filename = 'viscosityConstModel.dat'
    inquire(file=filename, exist=file_exists)
	if (file_exists) then
        open(fnum, file='viscosityConstModel.dat', status='old', iostat=ierr)         
      	if (ierr .eq. int(0)) then           
            read(fnum,*,iostat=rerr) typeViscConstModel
        endif
        write(*,*) "typeViscConstModel =",typeViscConstModel
        !---------------------------------------------
        !Define type of Constitutive model
        !   1/Newtonian, nparam = 1
        !   2/Power-law model, nparam = 2
        !   3/Powell-Eyring model, nparam =3
        !   4/Cross model, nparam = 4
        !   5/Carreau model, nparam = 4
        !   6/Carreau-Yasuda model, nparam = 5
        !---------------------------------------------
        if (typeViscConstModel.eq.1) then
            nparamViscConstModel = 1
        elseif(nparamViscConstModel.eq.2) then
            nparamViscConstModel = 2
        elseif(nparamViscConstModel.eq.3) then
            nparamViscConstModel = 3
        elseif(nparamViscConstModel.eq.4) then
            nparamViscConstModel = 4
        elseif(nparamViscConstModel.eq.5) then
            nparamViscConstModel = 4
        elseif(nparamViscConstModel.eq.6) then
            nparamViscConstModel = 5
        else
            write(*,*) 'Constitutive model not defined '
            stop
        endif            

        allocate(paramViscConstModel(nparamViscConstModel))

        do i=1,nparamViscConstModel
             read(fnum,*) paramViscConstModel(i)
        enddo

	else  
        write(*,*) 'Error, Viscosity Constitutive Model. File ', filename, ' not found. Exiting.'
        stop
	endif



    close(fnum) 
endif    


end subroutine initialize_NNW


subroutine get_mu(rmu,gamma,npro)
implicit none
integer :: npro
real*8, intent(inout)  :: rmu(npro)
real*8, intent(inout)  :: gamma(npro,3)
real*8 :: mu, mu_inf, mu0, lambda, n, m, a

if (nparamViscConstModel.eq.1) then ! Newtonian 
    mu = paramViscConstModel(1)
    rmu = mu
endif



end subroutine get_mu

subroutine get_shear_rate (gradv_x1, gradv_x2, gradv_x3, npro, gamma)
implicit none
integer :: npro
real*8, intent(inout)  :: gamma(npro), gradv_x1(npro,3)
real*8, intent(inout)  :: gradv_x2(npro,3), gradv_x3(npro,3)
real*8 :: rod_12(npro), rod_23(npro), rod_13(npro)
real*8 :: rod_11(npro), rod_22(npro), rod_33(npro)

     ! Components of the rate of deformation tensor

rod_12 =  0.5d0*( gradv_x2(:,1) + gradv_x1(:,2) )  
rod_23 =  0.5d0*( gradv_x3(:,2) + gradv_x2(:,3) )  
rod_13 =  0.5d0*( gradv_x1(:,3) + gradv_x3(:,1) )  
      
rod_11 =  gradv_x1(:,1) 
rod_22 =  gradv_x2(:,2)  
rod_33 =  gradv_x3(:,3)  

gamma = sqrt(4.0d0*(rod_12**2.0d0 + rod_23**2.0d0 + rod_13**2.0d0) & 
           + 2.0d0*(rod_11**2.0d0 + rod_22**2.0d0 + rod_33**2.0d0))

end subroutine get_shear_rate

end module nnw