module lateral_grid
  use fourier_transform
  use lege_poly
  implicit none
  
  type, public :: T_lateralGrid
    type(T_legep), public :: lgp
    type(T_fft),   public :: fourtrans
    
    contains
    
    procedure :: init_sub       => init_harmonics_sub
    procedure :: deallocate_sub => deallocate_harmonics_sub
    
    procedure :: harmsy_sub
    procedure :: harman_sub
    
  end type T_lateralGrid
  
  interface
    module subroutine init_harmonics_sub(this, jmax)
      class(T_lateralGrid), intent(inout) :: this
      integer,              intent(in)    :: jmax
    end subroutine init_harmonics_sub
    
    module subroutine deallocate_harmonics_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_harmonics_sub
    
    module subroutine harmsy_sub(this, cin, grid)
      class(T_lateralGrid), intent(in)  :: this
      complex(kind=dbl),    intent(in)  :: cin(*)
      real(kind=dbl),       intent(out) :: grid(32,*)
    end subroutine harmsy_sub
    
    module subroutine harman_sub(this, grid, cout)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(inout) :: grid(*)
      complex(kind=dbl),    intent(out)   :: cout(*)
    end subroutine harman_sub
  end interface
  
end module lateral_grid
