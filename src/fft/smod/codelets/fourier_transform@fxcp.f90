submodule (fourier_transform) fxcp
  implicit none; contains
  
  module procedure fxcpy
    integer :: i1, i2
    
    !GCC$ unroll 8
    !DIR$ unroll (8)
    do i2 = 1, m
      !$omp simd
      do i1 = 1, ndbl
        arr_to(i1,i2) = arr_from(i1,i2)
      end do
    end do
    
  end procedure fxcpy
  
end submodule fxcp