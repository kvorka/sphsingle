submodule (physical_grid) alloc
  implicit none; contains
  
  module procedure alloc_grid_sub
    
    this%c_grid = malloc( default_alig, nth * nph * size_d )
    call c_f_pointer( this%c_grid, this%tp, [nth,nph] )
    
  end procedure alloc_grid_sub
  
end submodule alloc