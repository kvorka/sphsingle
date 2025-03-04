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
    
    module pure subroutine fxzm2a(howmany, k, l, x, t)
      integer,        intent(in)    :: k, l, howmany
      real(kind=dbl), intent(in)    :: t(0:*)
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:1,0:k-1)
    end subroutine fxzm2a
    
    module pure subroutine fxzm2b(howmany, l, x)
      integer,        intent(in)    :: l, howmany
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:1)
    end subroutine fxzm2b
    
    module pure subroutine fxzm3a(howmany, k, l, x, t)
      integer,        intent(in)    :: k, l, howmany
      real(kind=dbl), intent(in)    :: t(0:*)
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:2,0:k-1)
    end subroutine fxzm3a
    
    module pure subroutine fxzm3b(howmany, l, x)
      integer,        intent(in)    :: l, howmany
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:2)
    end subroutine fxzm3b
    
    module pure subroutine fxzm4a(howmany, k, l, x, t)
      integer,        intent(in)    :: k, l, howmany
      real(kind=dbl), intent(in)    :: t(0:*)
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:3,0:k-1)
    end subroutine fxzm4a
    
    module pure subroutine fxzm4b(howmany, l, x)
      integer,        intent(in)    :: l, howmany
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:3)
    end subroutine fxzm4b
    
    module pure subroutine fxzm5a(howmany, k, l, x, t)
      integer,        intent(in)    :: k, l, howmany
      real(kind=dbl), intent(in)    :: t(2,0:*)
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:4,0:k-1)
    end subroutine fxzm5a
    
    module pure subroutine fxzm5b(howmany, l, x)
      integer,        intent(in)    :: l, howmany
      real(kind=dbl), intent(inout) :: x(howmany,2,l,0:4)
    end subroutine fxzm5b
  end interface
  
end module fourier_transform