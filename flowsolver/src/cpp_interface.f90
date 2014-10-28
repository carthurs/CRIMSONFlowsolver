module cpp_interface
    interface
    		subroutine giveflowpointertocpp(surfaceIndex, flowPointer) bind(c,name="giveflowpointertocpp")
    			use iso_c_binding
    			integer(c_int) :: surfaceIndex
    			type(c_ptr) :: flowPointer
    		end subroutine giveflowpointertocpp
    end interface
    
    interface
    		subroutine givepressurepointertocpp(surfaceIndex, pressPointer) bind(c,name="givepressurepointertocpp")
    			use iso_c_binding
    			integer(c_int) :: surfaceIndex
    			type(c_ptr) :: pressPointer
    		end subroutine givepressurepointertocpp
    end interface

    interface
    		subroutine callCppComputeAllImplicitCoeff_solve(timestepNumber) bind(c,name="callCppComputeAllImplicitCoeff_solve")
    			use iso_c_binding
    			integer(c_int) :: timestepNumber
    		end subroutine callCppComputeAllImplicitCoeff_solve
    end interface

    interface
    		subroutine callCppComputeAllImplicitCoeff_update(timestepNumber) bind(c,name="callCppComputeAllImplicitCoeff_update")
    			use iso_c_binding
    			integer(c_int) :: timestepNumber
    		end subroutine callCppComputeAllImplicitCoeff_update
    end interface

    interface
    		subroutine callCppGetImplicitCoeff_rcr(implicitCoeffs_toBeFilled_ptr) bind(c,name="callCppGetImplicitCoeff_rcr")
    			use iso_c_binding
    			type(c_ptr) :: implicitCoeffs_toBeFilled_ptr
    		end subroutine callCppGetImplicitCoeff_rcr
    end interface


    interface
    		subroutine callCPPUpdateAllRCRS_Pressure_n1_withflow() bind(c,name="callCPPUpdateAllRCRS_Pressure_n1_withflow")
    			use iso_c_binding
    		end subroutine callCPPUpdateAllRCRS_Pressure_n1_withflow
    end interface

    interface
    		subroutine callCPPUpdateAllRCRS_setflow_n(numberOfRCRSurfaces, flows) bind(c,name="callCPPUpdateAllRCRS_setflow_n")
    			use iso_c_binding
    			integer(c_int) :: numberOfRCRSurfaces
    			type(c_ptr) :: flows
    		end subroutine callCPPUpdateAllRCRS_setflow_n
    end interface

    interface
    		subroutine callCPPUpdateAllRCRS_setflow_n1(numberOfRCRSurfaces, flows) bind(c,name="callCPPUpdateAllRCRS_setflow_n1")
    			use iso_c_binding
    			integer(c_int) :: numberOfRCRSurfaces
    			type(c_ptr) :: flows
    		end subroutine callCPPUpdateAllRCRS_setflow_n1
    end interface
    
end module cpp_interface