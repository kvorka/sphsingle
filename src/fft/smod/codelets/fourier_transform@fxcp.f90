submodule (fourier_transform) fxcp
  implicit none; contains
  
  module procedure fxcpy
    integer :: i1, i2
    
    do i2 = 1, m, 8
      !$omp simd
      do i1 = 1, ndbl
        arr_to(i1,i2  ) = arr_from(i1,i2  )
        arr_to(i1,i2+1) = arr_from(i1,i2+1)
        arr_to(i1,i2+2) = arr_from(i1,i2+2)
        arr_to(i1,i2+3) = arr_from(i1,i2+3)
        arr_to(i1,i2+4) = arr_from(i1,i2+4)
        arr_to(i1,i2+5) = arr_from(i1,i2+5)
        arr_to(i1,i2+6) = arr_from(i1,i2+6)
        arr_to(i1,i2+7) = arr_from(i1,i2+7)
      end do
    end do
    
  end procedure fxcpy
  
end submodule fxcp