submodule (fourier_transform) fxsc
  implicit none; contains
  
  module procedure fxrsc
    integer :: i1, i2
    
    do i2 = 1, m, 2
      !$omp simd
      do i1 = 1, ndbl
        arr(i1,i2  ) = fac * arr(i1,i2  )
        arr(i1,i2+1) = fac * arr(i1,i2+1)
      end do
    end do
    
  end procedure fxrsc
  
end submodule fxsc