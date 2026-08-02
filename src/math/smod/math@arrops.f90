submodule (math) arrops
  implicit none; contains
  
  module procedure zero_rarray_sub
    integer :: i1, i2
    
    !! Vectorize the inner loop and unroll the outer loop.
    !GCC$ unroll 16
    !DIR$ unroll (16)
    do i2 = 1, (n/4)*4, 4
      !$omp simd
      do i1 = 0, 3
        arr(i1+i2) = 0._dbl
      end do
    end do
    
  end procedure zero_rarray_sub
  
end submodule arrops