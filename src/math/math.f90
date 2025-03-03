module math
  use iso_fortran_env, only: real64, real128
  use iso_c_binding
  implicit none; public
  
  integer, parameter :: dbl  = real64   !double precision
  integer, parameter :: qbl  = real128  !quadruple precision
  
  interface
    module pure subroutine zero_rarray_sub(n, arr)
      integer,        intent(in)  :: n
      real(kind=dbl), intent(out) :: arr(n)
    end subroutine zero_rarray_sub
  end interface
  
  interface
    type(c_ptr) function malloc(alignement, n) bind(C, name='aligned_alloc')
      import         :: c_ptr
      integer, value :: alignement, n
    end function malloc
    
    subroutine free(ptr) bind(C, name="free")
      import             :: c_ptr
      type(c_ptr), value :: ptr
    end subroutine free
  end interface
  
end module math
