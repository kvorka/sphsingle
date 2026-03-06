module memloc
  use iso_fortran_env, only: qbl => real128, dbl => real64
  use iso_c_binding,   only: c_ptr, c_sizeof, c_f_pointer
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
  
end module memloc