submodule (fourier_transform) fx5b
  implicit none; contains
  
  module procedure fxzm5b_c
    integer        :: i, iv, iv1
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im, x4re, x4im
    
    do i = 1, l
      do iv = 1, howmany, 16
        !$omp simd
        do iv1 = 0, 15
          x0re = x(iv1+iv,1,i,1) - x(iv1+iv,1,i,4)
          x0im = x(iv1+iv,2,i,1) - x(iv1+iv,2,i,4)
          x1re = x(iv1+iv,1,i,1) + x(iv1+iv,1,i,4)
          x1im = x(iv1+iv,2,i,1) + x(iv1+iv,2,i,4)
          x3re = x(iv1+iv,1,i,2) - x(iv1+iv,1,i,3)
          x3im = x(iv1+iv,2,i,2) - x(iv1+iv,2,i,3)
          x4re = x(iv1+iv,1,i,2) + x(iv1+iv,1,i,3)
          x4im = x(iv1+iv,2,i,2) + x(iv1+iv,2,i,3)
          
          x2re = C53 * x3re + x0re
          x2im = C53 * x3im + x0im
          x3re = C53 * x0re - x3re
          x3im = C53 * x0im - x3im
          x0re =       x1re + x4re
          x0im =       x1im + x4im
          x1re =       x1re - x4re
          x1im =       x1im - x4im
          
          x4re =     x(iv1+iv,1,i,0) - C51 * x0re
          x4im =     x(iv1+iv,2,i,0) - C51 * x0im
          x1re =     x4re            - C52 * x1re
          x1im =     x4im            - C52 * x1im
          x4re = 2 * x4re            -       x1re
          x4im = 2 * x4im            -       x1im
          
          x(iv1+iv,1,i,0) =     x(iv1+iv,1,i,0) +       x0re
          x(iv1+iv,2,i,0) =     x(iv1+iv,2,i,0) +       x0im
          x(iv1+iv,1,i,3) =     x1re            - C54 * x3im
          x(iv1+iv,2,i,3) =     x1im            + C54 * x3re
          x(iv1+iv,1,i,2) = 2 * x1re            -       x(iv1+iv,1,i,3)
          x(iv1+iv,2,i,2) = 2 * x1im            -       x(iv1+iv,2,i,3)
          x(iv1+iv,1,i,4) =     x4re            - C54 * x2im
          x(iv1+iv,2,i,4) =     x4im            + C54 * x2re
          x(iv1+iv,1,i,1) = 2 * x4re            -       x(iv1+iv,1,i,4)
          x(iv1+iv,2,i,1) = 2 * x4im            -       x(iv1+iv,2,i,4)
        end do
      end do
    end do
    
  end procedure fxzm5b_c
  
end submodule fx5b