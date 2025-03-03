submodule (lege_poly) allocators
  implicit none; contains
  
  module procedure allocate_lgp_arr_sub
    integer :: n
    
    n = 4*this%nrma
    
    c_arr = malloc( 32, int( n * c_sizeof(0._dbl) ) )
    call c_f_pointer( c_arr, arr, [n] )
    call zero_rarray_sub( n, arr )
    
  end procedure allocate_lgp_arr_sub
  
end submodule allocators