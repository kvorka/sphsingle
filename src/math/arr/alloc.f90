submodule (math) alloc
  implicit none
  
  integer, parameter :: size_d = c_sizeof(0._dbl)
  
  contains
  
  module procedure alloc_aligned_sub
    
    c_arr = malloc( alig, n * size_d )
      call c_f_pointer( c_arr, f_arr, [n] )
    
    call zero_rarray_c( n, f_arr )
    
  end procedure alloc_aligned_sub
  
end submodule alloc