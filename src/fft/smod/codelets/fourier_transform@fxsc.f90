submodule (fourier_transform) fxsc
  implicit none; contains
  
  module procedure fxrsc
    integer :: i1, i2
    
    !GCC$ unroll 4
    !DIR$ unroll (4)
    do i2 = 1, m
      !$omp simd
      do i1 = 1, ndbl
        arr(i1,i2) = fac * arr(i1,i2)
      end do
    end do
    
  end procedure fxrsc
  
end submodule fxsc