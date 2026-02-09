submodule (physical_grid) alloc
  implicit none; contains
  
  module procedure alloc_grid_sub
    
    call alloc_aligned_sub( nth*nph, this%c_grid, this%grid )
    
    this%tp(1:nth,1:nph) => this%grid(1:nth*nph)
    
  end procedure alloc_grid_sub
  
end submodule alloc