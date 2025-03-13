submodule (fourier_transform) fx2a
  implicit none; contains
  
  module procedure fxzm2a_c
    integer        :: i, j, iv, iv1
    real(kind=dbl) :: x0re, x0im, x1re, x1im, t1re, t1im
    
    do j = 0, k-1
      t1re = t(  2*j)
      t1im = t(1+2*j)
      
      do i = 1, l
        do iv = 1, howmany, 16
          !$omp simd
          do iv1 = 0, 15
            x0re = x(iv1+iv,1,i,0,j)
            x0im = x(iv1+iv,2,i,0,j)
            
            x1re = x0re - t1re * x(iv1+iv,1,i,1,j)
            x1im = x0im - t1im * x(iv1+iv,1,i,1,j)
            
            x(iv1+iv,1,i,1,j) =     x1re + t1im * x(iv1+iv,2,i,1,j)
            x(iv1+iv,2,i,1,j) =     x1im - t1re * x(iv1+iv,2,i,1,j)
            x(iv1+iv,1,i,0,j) = 2 * x0re -        x(iv1+iv,1,i,1,j)
            x(iv1+iv,2,i,0,j) = 2 * x0im -        x(iv1+iv,2,i,1,j)
          end do
        end do
      end do
    end do
    
  end procedure fxzm2a_c
  
end submodule fx2a
