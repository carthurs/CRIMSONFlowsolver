module cpp_interface
    interface
    		subroutine giveflowpointertocpp(surfaceIndex, flowPointer) bind(c,name="giveflowpointertocpp")
    			use iso_c_binding
    			integer(c_int) :: surfaceIndex
    			type(c_ptr) :: flowPointer
    		end subroutine
    end interface
end module cpp_interface