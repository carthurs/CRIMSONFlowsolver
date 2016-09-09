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
            subroutine callCPPGiveBoundaryConditionsListsOfTheirAssociatedMeshNodes(ndsurf_nodeToBoundaryAssociationArray, lengthOfNodeToBoundaryAssociationArray) bind(c,name="callCPPGiveBoundaryConditionsListsOfTheirAssociatedMeshNodes")
                use iso_c_binding
                type(c_ptr) :: ndsurf_nodeToBoundaryAssociationArray
                integer(c_int) :: lengthOfNodeToBoundaryAssociationArray
            end subroutine callCPPGiveBoundaryConditionsListsOfTheirAssociatedMeshNodes
    end interface


    interface
            subroutine callCPPGetBinaryMaskToAdjustNodalBoundaryConditions(binaryMask, binaryMaskLength) bind(c,name="callCPPGetBinaryMaskToAdjustNodalBoundaryConditions")
                use iso_c_binding
                type(c_ptr) :: binaryMask
                integer(c_int) :: binaryMaskLength
            end subroutine callCPPGetBinaryMaskToAdjustNodalBoundaryConditions
    end interface

    interface
            subroutine callCPPGetNumberOfBoundaryConditionsWhichCurrentlyDisallowFlow(numBCsWhichDisallowFlow) bind(c,name="callCPPGetNumberOfBoundaryConditionsWhichCurrentlyDisallowFlow")
                use iso_c_binding
                integer(c_int) :: numBCsWhichDisallowFlow
            end subroutine callCPPGetNumberOfBoundaryConditionsWhichCurrentlyDisallowFlow
    end interface

    interface
            subroutine callCPPGetNumberOfNetlistsWhichCurrentlyAllowFlow(numBCsWhichAllowFlow) bind(c,name="callCPPGetNumberOfNetlistsWhichCurrentlyAllowFlow")
                use iso_c_binding
                integer(c_int) :: numBCsWhichAllowFlow
            end subroutine callCPPGetNumberOfNetlistsWhichCurrentlyAllowFlow
    end interface

    interface
        subroutine callCPPGetNumberOfCppManagedBoundaryConditions(numManagedBCs) bind(c,name="callCPPGetNumberOfCppManagedBoundaryConditions")
            use iso_c_binding
            integer(c_int) :: numManagedBCs
        end subroutine callCPPGetNumberOfCppManagedBoundaryConditions
    end interface

    interface
            subroutine callCPPDiscoverWhetherFlowPermittedAcrossSurface(queriedSurfaceIndex,flowIsPermitted) bind(c,name="callCPPDiscoverWhetherFlowPermittedAcrossSurface")
                use iso_c_binding
                integer(c_int) :: queriedSurfaceIndex
                integer(c_int) :: flowIsPermitted
            end subroutine callCPPDiscoverWhetherFlowPermittedAcrossSurface
    end interface

    interface
            subroutine callCPPHaveBoundaryConditionTypesChanged(boundaryConditionTypesHaveChanged) bind(c,name="callCPPHaveBoundaryConditionTypesChanged")
                use iso_c_binding
                integer(c_int) :: boundaryConditionTypesHaveChanged
            end subroutine callCPPHaveBoundaryConditionTypesChanged
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
            subroutine callCPPSetPressureFromFortran() bind(c,name="callCPPSetPressureFromFortran")
                use iso_c_binding
            end subroutine callCPPSetPressureFromFortran
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
    ! interface
    !         subroutine callCppSetSurfacePressure_controlledCoronary(coronarySurfacePressures) bind(c,name="callCppSetSurfacePressure_controlledCoronary")
    !             use iso_c_binding
    !             type(c_ptr) :: coronarySurfacePressures
    !         end subroutine callCppSetSurfacePressure_controlledCoronary
    ! end interface

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
            subroutine callCppfinalizeLPNAtEndOfTimestep_controlledCoronary() bind(c,name="callCppfinalizeLPNAtEndOfTimestep_controlledCoronary")
                use iso_c_binding
            end subroutine callCppfinalizeLPNAtEndOfTimestep_controlledCoronary
    end interface

    interface
            subroutine callCppfinalizeLPNAtEndOfTimestep_netlists() bind(c,name="callCppfinalizeLPNAtEndOfTimestep_netlists")
                use iso_c_binding
            end subroutine callCppfinalizeLPNAtEndOfTimestep_netlists
    end interface
    

    ! interface
    !         subroutine callCPPUpdateAllControlledCoronaryLPNs_Pressure_n1_withflow() bind(c,name="callCPPUpdateAllControlledCoronaryLPNs_Pressure_n1_withflow")
    !             use iso_c_binding
    !         end subroutine callCPPUpdateAllControlledCoronaryLPNs_Pressure_n1_withflow
    ! end interface
    
!   ============= Controlled Coronary Block End =============
!   ============= Netlist LPN Block Start =============

    interface
            subroutine callCPPInitialiseLPNAtStartOfTimestep_netlist() bind(c,name="callCPPInitialiseLPNAtStartOfTimestep_netlist")
                use iso_c_binding
            end subroutine callCPPInitialiseLPNAtStartOfTimestep_netlist
    end interface

    interface
            subroutine callCPPUpdateAllNetlistLPNs(timestepNumber) bind(c,name="callCPPUpdateAllNetlistLPNs")
                use iso_c_binding
                integer(c_int) :: timestepNumber
            end subroutine callCPPUpdateAllNetlistLPNs
    end interface

    interface
            subroutine callCPPGetImplicitCoeff_netlistLPNs(implicitCoeffs_toBeFilled_ptr) bind(c,name="callCPPGetImplicitCoeff_netlistLPNs")
                use iso_c_binding
                type(c_ptr) :: implicitCoeffs_toBeFilled_ptr
            end subroutine callCPPGetImplicitCoeff_netlistLPNs
    end interface

    interface
            subroutine callCPPWriteAllNetlistComponentFlowsAndNodalPressures() bind(c,name="callCPPWriteAllNetlistComponentFlowsAndNodalPressures")
                use iso_c_binding
            end subroutine callCPPWriteAllNetlistComponentFlowsAndNodalPressures
    end interface

    interface
            subroutine callCPPDebugPrintFlowPointerTarget_BCM() bind(c,name="callCPPDebugPrintFlowPointerTarget_BCM")
                use iso_c_binding
            end subroutine callCPPDebugPrintFlowPointerTarget_BCM
    end interface

    ! interface
    !         subroutine callCPPLoadAllNetlistComponentFlowsAndNodalPressures() bind(c,name="callCPPLoadAllNetlistComponentFlowsAndNodalPressures")
    !             use iso_c_binding
    !         end subroutine callCPPLoadAllNetlistComponentFlowsAndNodalPressures
    ! end interface

    ! interface
    !         subroutine callCppSetSurfacePressure_netlistLPNs(netlistSurfacePressures) bind(c,name="callCppSetSurfacePressure_netlistLPNs")
    !             use iso_c_binding
    !             type(c_ptr) :: netlistSurfacePressures
    !         end subroutine callCppSetSurfacePressure_netlistLPNs
    ! end interface

!   ============= Netlist LPN Block End ===================
!   ============= Control Systems Block Start =============

    interface
            subroutine callCPPUpdateBoundaryConditionControlSystems() bind(c,name="callCPPUpdateBoundaryConditionControlSystems")
                use iso_c_binding
            end subroutine callCPPUpdateBoundaryConditionControlSystems
    end interface

    interface
        subroutine callCPPResetStateUsingKalmanFilteredEstimate(flow, pressure, surfaceIndex, timestepNumber) bind(c,name="callCPPResetStateUsingKalmanFilteredEstimate")
            use iso_c_binding
            double precision :: flow
            double precision :: pressure
            integer(c_int) :: surfaceIndex
            integer(c_int) :: timestepNumber
        end subroutine callCPPResetStateUsingKalmanFilteredEstimate
    end interface
    

end module cpp_interface