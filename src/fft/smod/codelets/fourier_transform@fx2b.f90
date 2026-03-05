submodule (fourier_transform) fx2b
  implicit none; contains
  
  module procedure fxzm2b
    integer :: i1, i2, i3
    
    !$omp parallel do private (i1,i2)
    do i3 = 1, l/2
      do i2 = 1, m
        !$omp simd
        do i1 = 1, ndbl
          x(i1,i2,1,i3,1) =     x(i1,i2,1,i3,0) - x(i1,i2,1,i3,1)
          x(i1,i2,2,i3,1) =     x(i1,i2,2,i3,0) - x(i1,i2,2,i3,1)
          x(i1,i2,1,i3,0) = 2 * x(i1,i2,1,i3,0) - x(i1,i2,1,i3,1)
          x(i1,i2,2,i3,0) = 2 * x(i1,i2,2,i3,0) - x(i1,i2,2,i3,1)
        end do
      end do
    end do
    !$omp end parallel do
    
  end procedure fxzm2b
  
end submodule fx2b
