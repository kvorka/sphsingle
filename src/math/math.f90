module math
  use memloc
  implicit none
  
  real(kind=dbl), parameter :: zero   = 0._dbl          !! double zero
  integer,        parameter :: size_d = c_sizeof(zero)  !! size of C double
  
  interface
    module subroutine alloc_aligned_sub( n, c_arr, f_arr )
      integer,                 intent(in)  :: n
      type(c_ptr),             intent(out) :: c_arr
      real(kind=dbl), pointer, intent(out) :: f_arr(:)
    end subroutine alloc_aligned_sub
    
    module subroutine zero_rarray_sub(n, arr)
      integer,        intent(in)  :: n
      real(kind=dbl), intent(out) :: arr(n)
    end subroutine zero_rarray_sub
    
    module real(kind=qbl) function lege_fn(deg, x)
      integer,        intent(in) :: deg
      real(kind=qbl), intent(in) :: x
    end function lege_fn
  end interface
  
end module math
