module math
  use iso_fortran_env, only: real64, real128
  implicit none; public
  
  integer, parameter :: dbl  = real64   !double precision
  integer, parameter :: qbl  = real128  !quadruple precision
  integer, parameter :: step = 32       !number of doubles handled at once
  
  interface
    module pure subroutine zero_rarray_sub(n, arr)
      integer,        intent(in)  :: n
      real(kind=dbl), intent(out) :: arr(n)
    end subroutine zero_rarray_sub
  end interface
  
end module math
