submodule (math) alloc
  implicit none; contains
  
  module procedure alloc_aligned_sub
    integer :: i
    
    c_arr = malloc( alig, n * size_d )
      call c_f_pointer( c_arr, f_arr, [n] )
    
    do i = 1, n
      f_arr(i) = 0._dbl
    end do
    
  end procedure alloc_aligned_sub
  
end submodule alloc