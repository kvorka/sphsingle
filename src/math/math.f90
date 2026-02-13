module math
  use iso_fortran_env, only: dbl => real64, &
                           & qbl => real128
  use iso_c_binding, only: c_ptr,      &
                         & c_sizeof,   &
                         & c_f_pointer
  implicit none; public

#if defined ( mem64 )
  integer, parameter :: alig = 64  !! avx512 alignement
  integer, parameter :: ndbl = 8   !! number of doubles in one avx512 registry
#elif defined ( mem32 )
  integer, parameter :: alig = 32  !! avx2 alignement
  integer, parameter :: ndbl = 4   !! number of doubles in one avx registry
#else
  integer, parameter :: alig = 16  !! default alignement fallback to SSE
  integer, parameter :: ndbl = 2   !! number of doubles in one SSE registry
#endif
  
  integer, parameter :: size_d = c_sizeof(0._dbl)  !! size of double
  integer, parameter :: step   = 2 * ndbl          !! stepping through the latitudinal grid
  
  interface
    module subroutine alloc_aligned_sub( n, c_arr, f_arr )
      integer,                 intent(in)  :: n
      type(c_ptr),             intent(out) :: c_arr
      real(kind=dbl), pointer, intent(out) :: f_arr(:)
    end subroutine alloc_aligned_sub
    
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
