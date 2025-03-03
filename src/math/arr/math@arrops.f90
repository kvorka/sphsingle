submodule (math) arrops
  implicit none; contains
  
  module procedure zero_rarray_sub
    integer :: i
    
    do i = 1, n
      arr(i) = 0._dbl
    end do
    
  end procedure zero_rarray_sub
  
end submodule arrops