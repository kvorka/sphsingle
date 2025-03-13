submodule (fourier_transform) fx2b
  implicit none; contains
  
  module procedure fxzm2b_c
    integer :: i, iv, iv1
    
    do i = 1, l
      do iv = 1, howmany, 16
        !$omp simd
        do iv1 = 0, 15
          x(iv1+iv,1,i,1) =     x(iv1+iv,1,i,0) - x(iv1+iv,1,i,1)
          x(iv1+iv,2,i,1) =     x(iv1+iv,2,i,0) - x(iv1+iv,2,i,1)
          x(iv1+iv,1,i,0) = 2 * x(iv1+iv,1,i,0) - x(iv1+iv,1,i,1)
          x(iv1+iv,2,i,0) = 2 * x(iv1+iv,2,i,0) - x(iv1+iv,2,i,1)
        end do
      end do
    end do
    
  end procedure fxzm2b_c
  
end submodule fx2b
