module math
  use memloc
  implicit none
  
  interface
    module subroutine alloc_aligned_sub( n, c_arr, f_arr )
      integer,                 intent(in)  :: n
      type(c_ptr),             intent(out) :: c_arr
      real(kind=dbl), pointer, intent(out) :: f_arr(:)
    end subroutine alloc_aligned_sub
    
    module subroutine zero_rarray_sub(n, arr) bind(C, name="zero_rarray_c")
      integer, value, intent(in)  :: n
      real(kind=dbl), intent(out) :: arr(*)
    end subroutine zero_rarray_sub
    
    module real(kind=qbl) function lege_fn(deg, x)
      integer,        intent(in) :: deg
      real(kind=qbl), intent(in) :: x
    end function lege_fn
  end interface
  
end module math
