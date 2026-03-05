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
          x0re = x(i1,i2,1,i3,1) - x(i1,i2,1,i3,4)
          x0im = x(i1,i2,2,i3,1) - x(i1,i2,2,i3,4)
          x1re = x(i1,i2,1,i3,1) + x(i1,i2,1,i3,4)
          x1im = x(i1,i2,2,i3,1) + x(i1,i2,2,i3,4)
          x3re = x(i1,i2,1,i3,2) - x(i1,i2,1,i3,3)
          x3im = x(i1,i2,2,i3,2) - x(i1,i2,2,i3,3)
          x4re = x(i1,i2,1,i3,2) + x(i1,i2,1,i3,3)
          x4im = x(i1,i2,2,i3,2) + x(i1,i2,2,i3,3)
          
          x2re =       x0re + C53 * x3re
          x2im =       x0im + C53 * x3im
          x3re = C53 * x0re -       x3re
          x3im = C53 * x0im -       x3im
          x0re =       x1re +       x4re
          x0im =       x1im +       x4im
          x1re =       x1re -       x4re
          x1im =       x1im -       x4im
          
          x4re =     x(i1,i2,1,i3,0) - C51 * x0re
          x4im =     x(i1,i2,2,i3,0) - C51 * x0im
          x1re =     x4re            - C52 * x1re
          x1im =     x4im            - C52 * x1im
          x4re = 2 * x4re            -       x1re
          x4im = 2 * x4im            -       x1im
          
          x(i1,i2,1,i3,0) =     x(i1,i2,1,i3,0) +       x0re
          x(i1,i2,2,i3,0) =     x(i1,i2,2,i3,0) +       x0im
          x(i1,i2,1,i3,3) =     x1re            - C54 * x3im
          x(i1,i2,2,i3,3) =     x1im            + C54 * x3re
          x(i1,i2,1,i3,2) = 2 * x1re            -       x(i1,i2,1,i3,3)
          x(i1,i2,2,i3,2) = 2 * x1im            -       x(i1,i2,2,i3,3)
          x(i1,i2,1,i3,4) =     x4re            - C54 * x2im
          x(i1,i2,2,i3,4) =     x4im            + C54 * x2re
          x(i1,i2,1,i3,1) = 2 * x4re            -       x(i1,i2,1,i3,4)
          x(i1,i2,2,i3,1) = 2 * x4im            -       x(i1,i2,2,i3,4)
        end do
      end do
    end do
    !$omp end parallel do
    
  end procedure fxzm5b
  
end submodule fx5b
