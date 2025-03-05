module math
  use iso_fortran_env, only: real64, real128
  use iso_c_binding
  implicit none; public
  
  integer, parameter :: dbl          = real64   !double precision
  integer, parameter :: qbl          = real128  !quadruple precision
  integer, parameter :: default_alig = 32       !memory alignement: SSE (16), AVX (32), AVX512 (64)
  
  interface
    type(c_ptr) function malloc(alignement, n) bind(C, name='aligned_alloc')
      import         :: c_ptr
      integer, value :: alignement, n
    end function malloc
    
    module subroutine alloc_aligned_sub( alig, n, c_arr, f_arr )
      integer,                 intent(in)  :: alig, n
      type(c_ptr),             intent(out) :: c_arr
      real(kind=dbl), pointer, intent(out) :: f_arr(:)
    end subroutine alloc_aligned_sub
    
    subroutine free(ptr) bind(C, name="free")
      import             :: c_ptr
      type(c_ptr), value :: ptr
    end subroutine free
  end interface
  
#ifdef omp
  interface
    module pure subroutine zero_rarray_c(n, arr)
      integer,        intent(in)  :: n
      real(kind=dbl), intent(out) :: arr(n)
    end subroutine zero_rarray_c
  end interface
#else
  interface
    module pure subroutine zero_rarray_c(n, arr) bind(C, name="zero_rarray_c")
      integer, value, intent(in)  :: n
      real(kind=dbl), intent(out) :: arr(*)
    end subroutine zero_rarray_c
  end interface
#endif
  
end module math
