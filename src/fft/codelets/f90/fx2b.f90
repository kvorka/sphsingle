submodule (fourier_transform) fx2b
  implicit none; contains
  
  module procedure fxzm2b_c
    integer :: i, iv
    
    do i = 1, l
      do iv = 1, howmany
        x(iv,1,i,1) =     x(iv,1,i,0) - x(iv,1,i,1)
        x(iv,2,i,1) =     x(iv,2,i,0) - x(iv,2,i,1)
        x(iv,1,i,0) = 2 * x(iv,1,i,0) - x(iv,1,i,1)
        x(iv,2,i,0) = 2 * x(iv,2,i,0) - x(iv,2,i,1)
      end do
    end do
    
  end procedure fxzm2b_c
  
end submodule fx2b
