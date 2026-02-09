submodule (math) arrops
  implicit none; contains
  
  module procedure zero_rarray_c
    integer :: i2, i1
    
    do i2 = 1, (n/32)*32, 32
      !$omp simd
      do i1 = 0, 31
        arr(i1+i2) = 0._dbl
      end do
    end do
    
    do i2 = (n/32)*32+1, n, 2
      !$omp simd
      do i1 = 0, 1
        arr(i1+i2) = 0._dbl
      end do
    end do
    
  end procedure zero_rarray_c
  
end submodule arrops