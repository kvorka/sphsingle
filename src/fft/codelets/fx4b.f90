submodule (fourier_transform) fx4b
  implicit none; contains
  
  module procedure fxzm4b_c
    integer        :: i, iv, iv1
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im
    
    do i = 1, l
      do iv = 1, howmany, 16
        !$omp simd
        do iv1 = 0, 15
          x2re = x(iv1+iv,1,i,0) - x(iv1+iv,1,i,2)
          x2im = x(iv1+iv,2,i,0) - x(iv1+iv,2,i,2)
          x0re = x(iv1+iv,1,i,0) + x(iv1+iv,1,i,2)
          x0im = x(iv1+iv,2,i,0) + x(iv1+iv,2,i,2)
          x3re = x(iv1+iv,1,i,1) - x(iv1+iv,1,i,3)
          x3im = x(iv1+iv,2,i,1) - x(iv1+iv,2,i,3)
          x1re = x(iv1+iv,1,i,1) + x(iv1+iv,1,i,3)
          x1im = x(iv1+iv,2,i,1) + x(iv1+iv,2,i,3)
          
          x(iv1+iv,1,i,2) = x0re - x1re
          x(iv1+iv,2,i,2) = x0im - x1im
          x(iv1+iv,1,i,0) = x0re + x1re
          x(iv1+iv,2,i,0) = x0im + x1im       
          x(iv1+iv,1,i,1) = x2re - x3im
          x(iv1+iv,2,i,1) = x2im + x3re
          x(iv1+iv,1,i,3) = x2re + x3im
          x(iv1+iv,2,i,3) = x2im - x3re
        end do
      end do
    end do
    
  end procedure fxzm4b_c
  
end submodule fx4b