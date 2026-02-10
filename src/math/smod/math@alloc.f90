submodule (math) alloc
  implicit none; contains
  
  module procedure alloc_aligned_sub
    
    c_arr = malloc( alig, n * size_d )
      call c_f_pointer( c_arr, f_arr, [n] )
    
  end procedure alloc_aligned_sub
  
end submodule alloc