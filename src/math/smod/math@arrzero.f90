submodule (math) arrzero
  implicit none; contains
  
  module procedure zero_rarray_sub
    integer :: i2, i1, n32
    
    !! Hand-unroll the loop by a factor of 32.
    n32 = (n/32)*32
    
    !! Main loop unrolled by 32. Vectorize the inner loop.
    do i2 = 1, n32, 32
      !$omp simd
      do i1 = 0, 31
        arr(i1+i2) = 0._dbl
      end do
    end do
    
    !! Safe to assume by the code structure, that the remainder
    !! is divisible by 4. Again, vectorize the inner loop.
    do i2 = n32+1, n, 4
      !$omp simd
      do i1 = 0, 3
        arr(i1+i2) = 0._dbl
      end do
    end do
    
  end procedure zero_rarray_sub
  
end submodule arrzero