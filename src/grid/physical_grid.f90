module physical_grid
  use math
  implicit none
  
  type, public :: T_grid
    type(c_ptr)                         :: c_grid
    real(kind=dbl), pointer, contiguous :: tp(:,:)
    
    contains
    
    procedure :: alloc_sub => alloc_grid_sub
    procedure :: free_sub  => free_grid_sub
    
  end type T_grid
  
  interface
    module subroutine alloc_grid_sub(this, nth, nph)
      class(T_grid), intent(inout) :: this
      integer,       intent(in)    :: nth, nph
    end subroutine alloc_grid_sub
    
    module subroutine free_grid_sub(this)
      class(T_grid), intent(inout) :: this
    end subroutine free_grid_sub
  end interface
  
end module physical_grid