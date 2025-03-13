submodule (math) alloc
  implicit none; contains
  
  module procedure alloc_aligned_sub
    
    c_arr = malloc( alig, n * size_d )
      call c_f_pointer( c_arr, f_arr, [n] )
    
    call zero_rarray_c( n, f_arr )
    
  end procedure alloc_aligned_sub
  
  module procedure alloc_aligned_2d_sub
    
    c_arr = malloc( alig, n2 * n1 * size_d )
      call c_f_pointer( c_arr, f_arr, [n1,n2] )
    
    call zero_rarray_c( n1 * n2, f_arr )
    
  end procedure alloc_aligned_2d_sub
  
end submodule alloc