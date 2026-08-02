submodule (fourier_transform) fx3b
  implicit none; contains
  
  module procedure fxzm3b
    integer        :: i1, i2, i3
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im
    
    !$omp parallel do private (i1,i2,x0re,x0im,x1re,x1im,x2re,x2im)
    do i3 = 1, l/3
      !GCC$ unroll 4
      !DIR$ unroll (4)
      do i2 = 1, m
        !$omp simd
        do i1 = 1, ndbl
          x1re = x(i1,i2,1,i3,1) -       x(i1,i2,1,i3,2)
          x1im = x(i1,i2,2,i3,1) -       x(i1,i2,2,i3,2)
          x0re = x(i1,i2,1,i3,1) +       x(i1,i2,1,i3,2)
          x0im = x(i1,i2,2,i3,1) +       x(i1,i2,2,i3,2)
          x2re = x(i1,i2,1,i3,0) + C31 * x0re
          x2im = x(i1,i2,2,i3,0) + C31 * x0im
          
          x(i1,i2,1,i3,0) =     x0re +       x(i1,i2,1,i3,0)
          x(i1,i2,2,i3,0) =     x0im +       x(i1,i2,2,i3,0)
          x(i1,i2,1,i3,2) =     x2re + C32 * x1im
          x(i1,i2,2,i3,2) =     x2im - C32 * x1re
          x(i1,i2,1,i3,1) = 2 * x2re -       x(i1,i2,1,i3,2)
          x(i1,i2,2,i3,1) = 2 * x2im -       x(i1,i2,2,i3,2)
        end do
      end do
    end do
    !$omp end parallel do
    
  end procedure fxzm3b
  
end submodule fx3b