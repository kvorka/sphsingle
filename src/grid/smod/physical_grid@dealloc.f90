submodule (physical_grid) dealloc
  implicit none; contains
  
  module procedure free_grid_sub
    
    this%tp => null()
    call free( this%c_grid )
    
  end procedure free_grid_sub
  
end submodule dealloc