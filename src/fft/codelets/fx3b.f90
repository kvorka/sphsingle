submodule (fourier_transform) fx3b
  implicit none; contains
  
  module procedure fxzm3b_c
    integer        :: i, iv, iv1
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im
    
    do i = 1, l
      do iv = 1, howmany, 16
        !$omp simd
        do iv1 = 0, 15
          x0re = x(iv1+iv,1,i,1) + x(iv1+iv,1,i,2)
          x0im = x(iv1+iv,2,i,1) + x(iv1+iv,2,i,2)
          x1re = x(iv1+iv,1,i,1) - x(iv1+iv,1,i,2)
          x1im = x(iv1+iv,2,i,1) - x(iv1+iv,2,i,2)
          x2re = x(iv1+iv,1,i,0) + C31 * x0re
          x2im = x(iv1+iv,2,i,0) + C31 * x0im
          
          x1re = C32 * x1re
          x1im = C32 * x1im
          
          x(iv1+iv,1,i,0) = x0re + x(iv1+iv,1,i,0)
          x(iv1+iv,2,i,0) = x0im + x(iv1+iv,2,i,0)
          x(iv1+iv,1,i,1) = x2re - x1im
          x(iv1+iv,2,i,1) = x2im + x1re
          x(iv1+iv,1,i,2) = x2re + x1im
          x(iv1+iv,2,i,2) = x2im - x1re
        end do
      end do
    end do
    
  end procedure fxzm3b_c
  
end submodule fx3b
