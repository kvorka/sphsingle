submodule (physical_grid) alloc
  implicit none; contains
  
  module procedure alloc_grid_sub
    
    call alloc_aligned_2d_sub( nth, nph, this%c_grid, this%tp )
    
  end procedure alloc_grid_sub
  
end submodule alloc