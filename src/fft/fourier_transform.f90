module fourier_transform
  !! Original work: fxpack (ISPACK FORTRAN SUBROUTINE LIBRARY FOR SCIENTIFIC COMPUTING, Keiichi Ishioka)
  use math
  implicit none
  
  type, public :: T_fft
    integer                     :: n
    integer,        allocatable :: it(:)
    real(kind=dbl), allocatable :: t(:)
    
    contains
    
    procedure, public, pass :: init_sub => fft_init_sub
    procedure, public, pass :: fft_bwd_sub
    procedure, public, pass :: fft_fwd_sub
    procedure, public, pass :: deallocate_sub => fft_deallocate_sub
    
  end type T_fft
  
  real(kind=dbl), parameter :: pi = acos(-1._dbl)
  
  interface
    module subroutine fft_init_sub(this, n)
      class(T_fft), intent(inout) :: this
      integer,      intent(in)    :: n
    end subroutine fft_init_sub
    
    module subroutine fft_deallocate_sub(this)
      class(T_fft), intent(inout) :: this
    end subroutine fft_deallocate_sub
    
    module subroutine fft_fwd_sub(this, m, x)
      class(T_fft),      intent(in)    :: this
      integer,           intent(in)    :: m
      real(kind=dbl),    intent(inout) :: x(4*m*ndbl,2,0:this%n/2-1)
    end subroutine fft_fwd_sub
    
    module subroutine fft_bwd_sub(this, m, x)
      class(T_fft),   intent(in)    :: this
      integer,        intent(in)    :: m
      real(kind=dbl), intent(inout) :: x(4*m*ndbl,2,0:this%n/2-1)
    end subroutine fft_bwd_sub
    
    module subroutine fxzini(n, it, t)
      integer,        intent(in)  :: n
      integer,        intent(out) :: it(n)
      real(kind=dbl), intent(out) :: t(2,0:n-1)
    end subroutine fxzini
    
    module subroutine fxrsc(m, fac, arr) bind(C, name="fxrsc_c")
      integer,        value, intent(in)    :: m
      real(kind=dbl), value, intent(in)    :: fac
      real(kind=dbl),        intent(inout) :: arr(*)
    end subroutine fxrsc
    
    module subroutine fxztal(n, it, t, m, x) bind(C, name="fxztal_c")
      integer, value, intent(in)    :: n, m
      integer,        intent(in)    :: it(*)
      real(kind=dbl), intent(in)    :: t(*)
      real(kind=dbl), intent(inout) :: x(*)
    end subroutine fxztal
    
    module subroutine fxrc0(m, arr1, arr2) bind(C, name="fxrc0_c")
      integer, value, intent(in)    :: m
      real(kind=dbl), intent(inout) :: arr1(*), arr2(*)
    end subroutine fxrc0
    
    module subroutine fxc2r(m, t, x11, x12, x21, x22) bind(C, name="fxc2r_c")
      integer, value, intent(in)    :: m
      real(kind=dbl), intent(in)    :: t(*)
      real(kind=dbl), intent(inout) :: x11(*), x12(*), x21(*), x22(*)
    end subroutine fxc2r
    
    module subroutine fxr2c(m, t, x11, x12, x21, x22) bind(C, name="fxr2c_c")
      integer, value, intent(in)    :: m
      real(kind=dbl), intent(in)    :: t(*)
      real(kind=dbl), intent(inout) :: x11(*), x12(*), x21(*), x22(*)
    end subroutine fxr2c
  end interface
  
end module fourier_transform
