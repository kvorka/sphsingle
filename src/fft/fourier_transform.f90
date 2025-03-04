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
    
    procedure, public,  pass :: init_sub       => fft_init_sub
    procedure, public,  pass :: deallocate_sub => fft_deallocate_sub
    
    procedure, private, pass :: fxztal
    procedure, private, pass :: fxzshf
    
    procedure, public,  pass :: fft_r2c_sub
    procedure, public,  pass :: fft_c2r_sub
    
  end type T_fft
  
  integer,        parameter :: imm = -2e4
  real(kind=dbl), parameter :: pi  = acos(-1._dbl)
  real(kind=dbl), parameter :: C31 = -0.5_dbl
  real(kind=dbl), parameter :: C32 = +0.86602540378443864676_dbl
  real(kind=dbl), parameter :: C51 = +0.25_dbl
  real(kind=dbl), parameter :: C52 = +0.5590169943749474241_dbl
  real(kind=dbl), parameter :: C53 = +0.6180339887498948482_dbl
  real(kind=dbl), parameter :: C54 = -0.9510565162951535721_dbl
  
  interface
    module pure subroutine fft_init_sub(this, n)
      class(T_fft), intent(inout) :: this
      integer,      intent(in)    :: n
    end subroutine fft_init_sub
    
    module pure subroutine fft_deallocate_sub(this)
      class(T_fft), intent(inout) :: this
    end subroutine fft_deallocate_sub
    
     module subroutine fft_r2c_sub(this, howmany, x)
      class(T_fft),   intent(in)    :: this
      integer,        intent(in)    :: howmany
      real(kind=dbl), intent(inout) :: x(howmany,2,0:this%n/2-1)
    end subroutine fft_r2c_sub
    
    module subroutine fft_c2r_sub(this, howmany, x)
      class(T_fft),   intent(in)    :: this
      integer,        intent(in)    :: howmany
      real(kind=dbl), intent(inout) :: x(howmany,2,0:this%n/2-1)
    end subroutine fft_c2r_sub
    
    module pure subroutine fxztal(this, howmany, x)
      class(T_fft),   intent(in)    :: this
      integer,        intent(in)    :: howmany
      real(kind=dbl), intent(inout) :: x(*)
    end subroutine fxztal
    
    module subroutine fxzshf(this, howmany, x)
      class(T_fft),   intent(in)    :: this
      integer,        intent(in)    :: howmany
      real(kind=dbl), intent(inout) :: x(2*howmany,0:this%n/2-1)
    end subroutine fxzshf
  end interface
  
  interface
    module pure function prime_adjustement_sub(n) result (nn)
      integer, intent(in) :: n
      integer             :: nn
    end function prime_adjustement_sub
    
    module pure subroutine fxzini(n, it, t)
      integer,        intent(in)  :: n
      integer,        intent(out) :: it(n)
      real(kind=dbl), intent(out) :: t(2,0:n-1)
    end subroutine fxzini
  end interface

#ifdef omp
  interface
    module pure subroutine fxzm2a_c(howmany, k, l, x, t)
      integer,        intent(in)    :: k, l, howmany
      real(kind=dbl), intent(in)    :: t(0:*)
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:1,0:k-1)
    end subroutine fxzm2a_c
    
    module pure subroutine fxzm2b_c(howmany, l, x)
      integer,        intent(in)    :: l, howmany
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:1)
    end subroutine fxzm2b_c
    
    module pure subroutine fxzm3a_c(howmany, k, l, x, t)
      integer,        intent(in)    :: k, l, howmany
      real(kind=dbl), intent(in)    :: t(0:*)
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:2,0:k-1)
    end subroutine fxzm3a_c
    
    module pure subroutine fxzm3b_c(howmany, l, x)
      integer,        intent(in)    :: l, howmany
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:2)
    end subroutine fxzm3b_c
    
    module pure subroutine fxzm4a_c(howmany, k, l, x, t)
      integer,        intent(in)    :: k, l, howmany
      real(kind=dbl), intent(in)    :: t(0:*)
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:3,0:k-1)
    end subroutine fxzm4a_c
    
    module pure subroutine fxzm4b_c(howmany, l, x)
      integer,        intent(in)    :: l, howmany
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:3)
    end subroutine fxzm4b_c
    
    module pure subroutine fxzm5a_c(howmany, k, l, x, t)
      integer,        intent(in)    :: k, l, howmany
      real(kind=dbl), intent(in)    :: t(2,0:*)
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:4,0:k-1)
    end subroutine fxzm5a_c
    
    module pure subroutine fxzm5b_c(howmany, l, x)
      integer,        intent(in)    :: l, howmany
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:4)
    end subroutine fxzm5b_c
  end interface
#else
  interface
    module pure subroutine fxzm2a_c(howmany, k, l, x, t)
      integer,        intent(in)    :: k, l, howmany
      real(kind=dbl), intent(in)    :: t(0:*)
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:1,0:k-1)
    end subroutine fxzm2a_c
    
    module pure subroutine fxzm2b_c(howmany, l, x)
      integer,        intent(in)    :: l, howmany
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:1)
    end subroutine fxzm2b_c
    
    module pure subroutine fxzm3a_c(howmany, k, l, x, t)
      integer,        intent(in)    :: k, l, howmany
      real(kind=dbl), intent(in)    :: t(0:*)
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:2,0:k-1)
    end subroutine fxzm3a_c
    
    module pure subroutine fxzm3b_c(howmany, l, x)
      integer,        intent(in)    :: l, howmany
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:2)
    end subroutine fxzm3b_c
    
    module pure subroutine fxzm4a_c(howmany, k, l, x, t) bind(C, name="fxzm4a_c")
      integer, value, intent(in)    :: k, l, howmany
      real(kind=dbl), intent(in)    :: t(*)
      real(kind=dbl), intent(inout) :: x(*)
    end subroutine fxzm4a_c
    
    module pure subroutine fxzm4b_c(howmany, l, x)
      integer,        intent(in)    :: l, howmany
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:3)
    end subroutine fxzm4b_c
    
    module pure subroutine fxzm5a_c(howmany, k, l, x, t)
      integer,        intent(in)    :: k, l, howmany
      real(kind=dbl), intent(in)    :: t(2,0:*)
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:4,0:k-1)
    end subroutine fxzm5a_c
    
    module pure subroutine fxzm5b_c(howmany, l, x)
      integer,        intent(in)    :: l, howmany
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:4)
    end subroutine fxzm5b_c
  end interface
#endif

end module fourier_transform