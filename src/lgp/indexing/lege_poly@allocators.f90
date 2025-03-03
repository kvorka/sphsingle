submodule (lege_poly) allocators
  implicit none; contains
  
  module procedure allocate_lgp_arr_sub
    integer :: n
    
    n = 4*this%nrma
    
    allocate( arr(n) )
      call zero_rarray_sub( n, arr )
    
  end procedure allocate_lgp_arr_sub
  
end submodule allocators