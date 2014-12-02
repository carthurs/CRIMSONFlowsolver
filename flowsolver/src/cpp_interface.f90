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


!   =============== Numerical RCR Block ===============
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
    		subroutine callCPPUpdateAllRCRS_setflow_n(flows) bind(c,name="callCPPUpdateAllRCRS_setflow_n")
    			use iso_c_binding
    			type(c_ptr) :: flows
    		end subroutine callCPPUpdateAllRCRS_setflow_n
    end interface

    interface
    		subroutine callCPPUpdateAllRCRS_setflow_n1(flows) bind(c,name="callCPPUpdateAllRCRS_setflow_n1")
    			use iso_c_binding
    			type(c_ptr) :: flows
    		end subroutine callCPPUpdateAllRCRS_setflow_n1
    end interface
!   ============= Numerical RCR Block End =============

    interface
            subroutine callCPPRecordPressuresAndFlowsInHistoryArrays() bind(c,name="callCPPRecordPressuresAndFlowsInHistoryArrays")
                use iso_c_binding
            end subroutine callCPPRecordPressuresAndFlowsInHistoryArrays
    end interface

    interface
            subroutine callCPPWritePHistAndQHistRCR() bind(c,name="callCPPWritePHistAndQHistRCR")
                use iso_c_binding
            end subroutine callCPPWritePHistAndQHistRCR
    end interface

!   =============== Controlled Coronary Block ===============
    interface
            subroutine callCppSetSurfacePressure_controlledCoronary(coronarySurfacePressures) bind(c,name="callCppSetSurfacePressure_controlledCoronary")
                use iso_c_binding
                type(c_ptr) :: coronarySurfacePressures
            end subroutine callCppSetSurfacePressure_controlledCoronary
    end interface

    interface
            subroutine callCppGetImplicitCoeff_controlledCoronary(implicitCoeffs_toBeFilled_ptr) bind(c,name="callCppGetImplicitCoeff_controlledCoronary")
                use iso_c_binding
                type(c_ptr) :: implicitCoeffs_toBeFilled_ptr
            end subroutine callCppGetImplicitCoeff_controlledCoronary
    end interface

    interface
            subroutine callCppUpdateAllControlledCoronaryLPNs() bind(c,name="callCppUpdateAllControlledCoronaryLPNs")
                use iso_c_binding
            end subroutine callCppUpdateAllControlledCoronaryLPNs
    end interface

    interface
            subroutine callCPPUpdateAllControlledCoronaryLPNs_Pressure_n1_withflow() bind(c,name="callCPPUpdateAllControlledCoronaryLPNs_Pressure_n1_withflow")
                use iso_c_binding
            end subroutine callCPPUpdateAllControlledCoronaryLPNs_Pressure_n1_withflow
    end interface
    
!   ============= Controlled Coronary Block End =============
    
end module cpp_interface