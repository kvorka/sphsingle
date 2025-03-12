submodule (lateral_grid) alloc
  implicit none; contains
  
  module procedure alloc_grid_sub
    
    c_grid = malloc( default_alig, 2 * this%lgp%n * this%fourtrans%n * size_d )
    call c_f_pointer( c_grid, f_grid, [2*this%lgp%n,this%fourtrans%n] )
    
  end procedure alloc_grid_sub
  
end submodule alloc