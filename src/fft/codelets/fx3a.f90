submodule (fourier_transform) fx3a
  implicit none; contains
  
  module procedure fxzm3a_c
    integer        :: i, j, iv, iv1
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im, t1re, t1im, t2re, t2im
    
    do j = 0, k-1
      t1re = t(  4*j)
      t1im = t(1+4*j)
      t2re = t(2+4*j)
      t2im = t(3+4*j)
      
      do i = 1, l
        do iv = 1, howmany, 16
          !$omp simd
          do iv1 = 0, 15
            x0re =        t1re * x(iv1+iv,1,i,1,j) - t1im * x(iv1+iv,2,i,1,j)
            x0im =        t1re * x(iv1+iv,2,i,1,j) + t1im * x(iv1+iv,1,i,1,j)
            x1re = x0re - t2re * x(iv1+iv,1,i,2,j) + t2im * x(iv1+iv,2,i,2,j)
            x1im = x0im - t2re * x(iv1+iv,2,i,2,j) - t2im * x(iv1+iv,1,i,2,j)
            
            x0re = 2 * x0re              -       x1re
            x0im = 2 * x0im              -       x1im
            x2re =     x(iv1+iv,1,i,0,j) + C31 * x0re
            x2im =     x(iv1+iv,2,i,0,j) + C31 * x0im
            
            x(iv1+iv,1,i,0,j) =     x0re +       x(iv1+iv,1,i,0,j)
            x(iv1+iv,2,i,0,j) =     x0im +       x(iv1+iv,2,i,0,j)
            x(iv1+iv,1,i,2,j) =     x2re + C32 * x1im
            x(iv1+iv,2,i,2,j) =     x2im - C32 * x1re
            x(iv1+iv,1,i,1,j) = 2 * x2re -       x(iv1+iv,1,i,2,j)
            x(iv1+iv,2,i,1,j) = 2 * x2im -       x(iv1+iv,2,i,2,j)
          end do
        end do
      end do
    end do
    
  end procedure fxzm3a_c
  
end submodule fx3a
