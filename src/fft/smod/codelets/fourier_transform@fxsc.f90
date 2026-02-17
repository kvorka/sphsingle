submodule (fourier_transform) fxsc
  implicit none; contains
  
  module procedure fxrsc
    integer :: i1, i2
    
    do i2 = 1, m, 4
      !$omp simd
      do i1 = 1, ndbl
        arr(i1,i2  ) = fac * arr(i1,i2  )
        arr(i1,i2+1) = fac * arr(i1,i2+1)
        arr(i1,i2+2) = fac * arr(i1,i2+2)
        arr(i1,i2+3) = fac * arr(i1,i2+3)
      end do
    end do
    
  end procedure fxrsc
  
end submodule fxsc