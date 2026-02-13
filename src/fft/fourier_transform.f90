module fourier_transform
  !Author of the original code: Keiichi Ishioka
  !Original work: fxpack (ISPACK FORTRAN SUBROUTINE LIBRARY FOR SCIENTIFIC COMPUTING)
  use math
  implicit none
  
  type, public :: T_fft
    integer                     :: n
    integer,        allocatable :: it(:)
    real(kind=dbl), allocatable :: t(:)
    
    contains
    
    procedure, public, pass :: init_sub       => fft_init_sub
    procedure, public, pass :: deallocate_sub => fft_deallocate_sub
    procedure, public, pass :: fft_r2c_sub, fft_c2r_sub
    
  end type T_fft
  
  integer,        parameter :: imm = -2e4
  real(kind=dbl), parameter :: pi  = acos(-1._dbl)
  
  interface
    module subroutine fft_init_sub(this, n)
      class(T_fft), intent(inout) :: this
      integer,      intent(in)    :: n
    end subroutine fft_init_sub
    
    module subroutine fft_deallocate_sub(this)
      class(T_fft), intent(inout) :: this
    end subroutine fft_deallocate_sub
    
     module subroutine fft_r2c_sub(this, m, x)
      class(T_fft),      intent(in)    :: this
      integer,           intent(in)    :: m
      real(kind=dbl),    intent(inout) :: x(m*ndbl,2,0:this%n/2-1)
    end subroutine fft_r2c_sub
    
    module subroutine fft_c2r_sub(this, m, x)
      class(T_fft),   intent(in)    :: this
      integer,        intent(in)    :: m
      real(kind=dbl), intent(inout) :: x(m*ndbl,2,0:this%n/2-1)
    end subroutine fft_c2r_sub
  end interface
  
  interface
    module function prime_adjustement_sub(n) result (nn)
      integer, intent(in) :: n
      integer             :: nn
    end function prime_adjustement_sub
    
    module subroutine fxzini(n, it, t)
      integer,        intent(in)  :: n
      integer,        intent(out) :: it(n)
      real(kind=dbl), intent(out) :: t(2,0:n-1)
    end subroutine fxzini
    
    module subroutine fxzshf(n, it, m, x)
      integer,        intent(in)    :: n, m, it(*)
      real(kind=dbl), intent(inout) :: x(m*ndbl,0:n/2-1)
    end subroutine fxzshf
    
    module subroutine fxztal(n, it, t, m, x)
      integer,        intent(in)    :: n, m, it(2)
      real(kind=dbl), intent(in)    :: t(0:*)
      real(kind=dbl), intent(inout) :: x(*)
    end subroutine fxztal
    
    module subroutine fxzm2a(m, k, l, x, t)
      integer,        intent(in)    :: m, k, l
      real(kind=dbl), intent(in)    :: t(2,0:*)
      real(kind=dbl), intent(inout) :: x(ndbl,m,2,l/2,0:1,0:k-1)
    end subroutine fxzm2a
    
    module subroutine fxzm2b(m, l, x)
      integer,        intent(in)    :: m, l
      real(kind=dbl), intent(inout) :: x(ndbl,m,2,l/2,0:1)
    end subroutine fxzm2b
    
    module subroutine fxzm3a(m, k, l, x, t)
      integer,        intent(in)    :: m, k, l
      real(kind=dbl), intent(in)    :: t(2,0:*)
      real(kind=dbl), intent(inout) :: x(ndbl,m,2,l/3,0:2,0:k-1)
    end subroutine fxzm3a
    
    module subroutine fxzm3b(m, l, x)
      integer,        intent(in)    :: m, l
      real(kind=dbl), intent(inout) :: x(ndbl,m,2,l/3,0:2)
    end subroutine fxzm3b
    
    module subroutine fxzm4a(m, k, l, x, t)
      integer,        intent(in)    :: m, k, l
      REAL(kind=dbl), intent(in)    :: t(2,0:*)
      real(kind=dbl), intent(inout) :: x(ndbl,m,2,l/4,0:3,0:k-1)
    end subroutine fxzm4a
    
    module subroutine fxzm4b(m, l, x)
      integer,        intent(in)    :: m, l
      real(kind=dbl), intent(inout) :: x(ndbl,m,2,l/4,0:3)
    end subroutine fxzm4b
    
    module subroutine fxzm5a(m, k, l, x, t)
      integer,        intent(in)    :: m, k, l
      real(kind=dbl), intent(in)    :: t(2,0:*)
      real(kind=dbl), intent(inout) :: x(ndbl,m,2,l/5,0:4,0:k-1)
    end subroutine fxzm5a
    
    module subroutine fxzm5b(m, l, x)
      integer,        intent(in)    :: m, l
      real(kind=dbl), intent(inout) :: x(ndbl,m,2,l/5,0:4)
    end subroutine fxzm5b
    
    module subroutine fxcpy(m, arr_from, arr_to)
      integer,        intent(in)  :: m
      real(kind=dbl), intent(in)  :: arr_from(ndbl,m)
      real(kind=dbl), intent(out) :: arr_to(ndbl,m)
    end subroutine fxcpy
    
    module subroutine fxaddsub(m, arr1, arr2)
      integer,        intent(in)    :: m
      real(kind=dbl), intent(inout) :: arr1(ndbl,m), arr2(ndbl,m)
    end subroutine fxaddsub
    
    module subroutine fxrsc(m, fac, arr)
      integer,        intent(in)    :: m
      integer,        intent(in)    :: fac
      real(kind=dbl), intent(inout) :: arr(ndbl,m)
    end subroutine fxrsc
    
    module subroutine fxc2r(m, t, x11, x12, x21, x22)
      integer,        intent(in)    :: m
      real(kind=dbl), intent(in)    :: t(2)
      real(kind=dbl), intent(inout) :: x11(ndbl,m), x12(ndbl,m), x21(ndbl,m), x22(ndbl,m)
    end subroutine fxc2r
    
    module subroutine fxr2c(m, t, x11, x12, x21, x22)
      integer,        intent(in)    :: m
      real(kind=dbl), intent(in)    :: t(2)
      real(kind=dbl), intent(inout) :: x11(ndbl,m), x12(ndbl,m), x21(ndbl,m), x22(ndbl,m)
    end subroutine fxr2c
  end interface
  
end module fourier_transform