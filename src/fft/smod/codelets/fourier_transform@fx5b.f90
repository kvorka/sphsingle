submodule (fourier_transform) fx5b
  implicit none; contains
  
  module procedure fxzm5b
    integer        :: i1, i2, i3
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im, x4re, x4im
    
    !$omp parallel do private (i1,i2,x0re,x0im,x1re,x1im,x2re,x2im,x3re,x3im,x4re,x4im)
    do i3 = 1, l/5
      do i2 = 1, m
        !$omp simd
        do i1 = 1, ndbl
          x0im = x(i1,i2,1,i3,1)
          x0re = x(i1,i2,1,i3,4)
          x2im = x(i1,i2,2,i3,1)
          x2re = x(i1,i2,2,i3,4)
          
          x1re = x0im + x0re
          x4re = x0im - x0re
          
          x0im = x(i1,i2,1,i3,2)
          x0re = x(i1,i2,1,i3,3)
          
          x1im = x2im + x2re
          x4im = x2im - x2re
          
          x2im = x(i1,i2,2,i3,2)
          x2re = x(i1,i2,2,i3,3)
          
          x3re = x0im - x0re
          x0re = x0im + x0re
          
          x3im = x2im - x2re
          x0im = x2im + x2re
          
          x2re = x4re + C53 * x3re
          x2im = x4im + C53 * x3im
          
          x3re = C53 * x4re - x3re
          x3im = C53 * x4im - x3im
          
          x4re = x1re + x0re
          x4im = x1im + x0im
          
          x1re = x1re - x0re
          x1im = x1im - x0im
          
          x0re = x(i1,i2,1,i3,0)
          x0im = x(i1,i2,2,i3,0)
          
          x(i1,i2,1,i3,0) = x0re + x4re
          x(i1,i2,2,i3,0) = x0im + x4im
          
          x0re = x0re - C51 * x4re
          x0im = x0im - C51 * x4im
          
          x1re = x0re - C52 * x1re
          x1im = x0im - C52 * x1im
          
          x0re = x0re + x0re
          x0im = x0im + x0im
          
          x0re = x0re - x1re
          x0im = x0im - x1im
          
          x3re = x1im + C54 * x3re
          x3im = x1re - C54 * x3im
          x2im = x0re - C54 * x2im
          x2re = x0im + C54 * x2re
          
          x(i1,i2,1,i3,4) = x2im
          x(i1,i2,2,i3,4) = x2re
          x(i1,i2,1,i3,3) = x3im
          x(i1,i2,2,i3,3) = x3re
          
          x1re = x1re + x1re
          x1im = x1im + x1im
          x0re = x0re + x0re
          x0im = x0im + x0im
          
          x1re = x1re - x3im
          x1im = x1im - x3re
          x0re = x0re - x2im
          x0im = x0im - x2re
          
          x(i1,i2,1,i3,2) = x1re
          x(i1,i2,2,i3,2) = x1im
          x(i1,i2,1,i3,1) = x0re
          x(i1,i2,2,i3,1) = x0im
        end do
      end do
    end do
    !$omp end parallel do
    
  end procedure fxzm5b
  
end submodule fx5b
