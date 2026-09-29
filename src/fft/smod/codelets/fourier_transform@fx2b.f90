submodule (fourier_transform) fx2b
  implicit none; contains
  
  module procedure fxzm2b
    integer        :: i1, i2, i3
    real(kind=dbl) :: x0re, x0im, x1re, x1im
    
    !$omp parallel do private (i1,i2,x0re,x0im,x1re,x1im)
    do i3 = 1, l/2
      do i2 = 1, m
        !$omp simd
        do i1 = 1, ndbl
          x0re = x(i1,i2,1,i3,0)
          x0im = x(i1,i2,2,i3,0)
          x1re = x(i1,i2,1,i3,1)
          x1im = x(i1,i2,2,i3,1)
          
          x1re = x0re - x1re
          x1im = x0im - x1im
          
          x0re = 2 * x0re - x1re
          x0im = 2 * x0im - x1im
          
          x(i1,i2,1,i3,1) = x1re
          x(i1,i2,2,i3,1) = x1im
          x(i1,i2,1,i3,0) = x0re
          x(i1,i2,2,i3,0) = x0im
        end do
      end do
    end do
    !$omp end parallel do
    
  end procedure fxzm2b
  
end submodule fx2b
