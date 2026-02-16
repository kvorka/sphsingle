submodule (fourier_transform) fx4b
  implicit none; contains
  
  module procedure fxzm4b
    integer        :: i1, i2, i3
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im
    
    do i3 = 1, l/4
      do i2 = 1, m
        !$omp simd
        do i1 = 1, ndbl
          x2re = x(i1,i2,1,i3,0) - x(i1,i2,1,i3,2)
          x2im = x(i1,i2,2,i3,0) - x(i1,i2,2,i3,2)
          x0re = x(i1,i2,1,i3,0) + x(i1,i2,1,i3,2)
          x0im = x(i1,i2,2,i3,0) + x(i1,i2,2,i3,2)
          x3re = x(i1,i2,1,i3,1) - x(i1,i2,1,i3,3)
          x3im = x(i1,i2,2,i3,1) - x(i1,i2,2,i3,3)
          x1re = x(i1,i2,1,i3,1) + x(i1,i2,1,i3,3)
          x1im = x(i1,i2,2,i3,1) + x(i1,i2,2,i3,3)
          
          x(i1,i2,1,i3,2) =     x0re - x1re
          x(i1,i2,2,i3,2) =     x0im - x1im
          x(i1,i2,1,i3,0) = 2 * x0re - x(i1,i2,1,i3,2)
          x(i1,i2,2,i3,0) = 2 * x0im - x(i1,i2,2,i3,2)       
          x(i1,i2,1,i3,1) =     x2re - x3im
          x(i1,i2,2,i3,1) =     x2im + x3re
          x(i1,i2,1,i3,3) = 2 * x2re - x(i1,i2,1,i3,1)
          x(i1,i2,2,i3,3) = 2 * x2im - x(i1,i2,2,i3,1)
        end do
      end do
    end do
    
  end procedure fxzm4b
  
end submodule fx4b