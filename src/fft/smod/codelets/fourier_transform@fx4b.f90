submodule (fourier_transform) fx4b
  implicit none; contains
  
  module procedure fxzm4b
    integer        :: i1, i2, i3
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im, x01, x02, x03, x04
    
    !$omp parallel do private (i1,i2,x0re,x0im,x1re,x1im,x2re,x2im,x3re,x3im,x01,x02,x03,x04)
    do i3 = 1, l/4
      do i2 = 1, m
        !$omp simd
        do i1 = 1, ndbl
          x01 = x(i1,i2,1,i3,0)
          x03 = x(i1,i2,1,i3,2)
          
          x0re = x01 + x03
          x2re = x01 - x03
          
          x02 = x(i1,i2,2,i3,0)
          x04 = x(i1,i2,2,i3,2)
          
          x0im = x02 + x04
          x2im = x02 - x04
          
          x01 = x(i1,i2,1,i3,1)
          x03 = x(i1,i2,1,i3,3)
          
          x1re = x01 + x03
          x3re = x01 - x03
          
          x02 = x(i1,i2,2,i3,1)
          x04 = x(i1,i2,2,i3,3)
          
          x1im = x02 + x04
          x3im = x02 - x04
          
          x1re = x0re - x1re
          x1im = x0im - x1im
          x01  = x0re + x0re
          x02  = x0im + x0im
          
          x(i1,i2,1,i3,2) = x1re
          x(i1,i2,2,i3,2) = x1im
          
          x0re = x01  - x1re
          x0im = x02  - x1im
          x3im = x2re - x3im
          x3re = x2im + x3re
          
          x01 = x2re + x2re
          x02 = x2im + x2im
          
          x(i1,i2,1,i3,0) = x0re
          x(i1,i2,2,i3,0) = x0im   
          x(i1,i2,1,i3,1) = x3im
          x(i1,i2,2,i3,1) = x3re
          
          x2re = x01 - x3im
          x2im = x02 - x3re
          
          x(i1,i2,1,i3,3) = x2re
          x(i1,i2,2,i3,3) = x2im
        end do
      end do
    end do
    !$omp end parallel do
    
  end procedure fxzm4b
  
end submodule fx4b