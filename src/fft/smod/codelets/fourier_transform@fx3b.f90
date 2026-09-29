submodule (fourier_transform) fx3b
  implicit none; contains
  
  module procedure fxzm3b
    integer        :: i1, i2, i3
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im, x01, x02, x03, x04
    
    !$omp parallel do private (i1,i2,x0re,x0im,x1re,x1im,x2re,x2im,x01,x02,x03,x04)
    do i3 = 1, l/3
      do i2 = 1, m
        !$omp simd
        do i1 = 1, ndbl
          x1re = x(i1,i2,1,i3,1)
          x2re = x(i1,i2,1,i3,2)
          
          x01 = x1re - x2re
          x03 = x1re + x2re
          
          x1im = x(i1,i2,2,i3,1)
          x2im = x(i1,i2,2,i3,2)
          
          x02 = x1im - x2im
          x04 = x1im + x2im
          
          x0re = x(i1,i2,1,i3,0)
          x0im = x(i1,i2,2,i3,0)
          
          x1re = x0re + C31 * x03
          x1im = x0im + C31 * x04
          
          x2re = x1re + C32 * x02
          x2im = x1im - C32 * x01
          
          x0re = x0re + x03
          x0im = x0im + x04
          x1re = x1re + x1re
          x1im = x1im + x1im
          
          x(i1,i2,1,i3,0) = x0re
          x(i1,i2,2,i3,0) = x0im
          x(i1,i2,1,i3,2) = x2re
          x(i1,i2,2,i3,2) = x2im
          
          x1re = x1re - x2re
          x1im = x1im - x2im
          
          x(i1,i2,1,i3,1) = x1re
          x(i1,i2,2,i3,1) = x1im
        end do
      end do
    end do
    !$omp end parallel do
    
  end procedure fxzm3b
  
end submodule fx3b