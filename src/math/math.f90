module math
  use iso_fortran_env, only: qbl => real128, dbl => real64
  use iso_c_binding,   only: c_ptr, c_f_pointer
  implicit none
  
#if defined (__AVX512F__)
  integer, parameter :: alig = 64  !! avx512 alignement
  integer, parameter :: ndbl = 8   !! number of doubles in one avx512 registry
#else
  integer, parameter :: alig = 32  !! avx2 alignement
  integer, parameter :: ndbl = 4   !! number of doubles in one avx registry
#endif 
  
  interface
    type(c_ptr) function malloc(alignement, n) bind(C, name='aligned_alloc')
      import         :: c_ptr
      integer, value :: alignement, n
    end function malloc
    
    subroutine free(ptr) bind(C, name="free")
      import             :: c_ptr
      type(c_ptr), value :: ptr
    end subroutine free
    
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
