submodule (math) arrops
  implicit none; contains
  
  module procedure zero_rarray_sub
    integer :: i1, i2
    
    !$omp simd
    do i2 = 1, n
      arr(i2) = 0._dbl
    end do
    
  end procedure zero_rarray_sub
  
end submodule arrops