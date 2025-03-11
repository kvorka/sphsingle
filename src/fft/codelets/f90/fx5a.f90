submodule (fourier_transform) fx5a
  implicit none; contains
  
  module procedure fxzm5a_c
    integer        :: i, j, ij, iv
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im, x4re, x4im, &
                    & t1re, t1im, t2re, t2im, t3re, t3im, t4re, t4im
    
    do j = 0, k-1
      ij = 4 * j
      
      t1re = t(1,ij  )
      t1im = t(2,ij  )
      t2re = t(1,ij+1)
      t2im = t(2,ij+1)
      t3re = t(1,ij+2)
      t3im = t(2,ij+2)
      t4re = t(1,ij+3)
      t4im = t(2,ij+3)
      
      do i = 1, l
        do iv = 1, howmany
          x1re =        t1re * x(iv,1,i,1,j) - t1im * x(iv,2,i,1,j)
          x1im =        t1re * x(iv,2,i,1,j) + t1im * x(iv,1,i,1,j)
          x2re =        t2re * x(iv,1,i,2,j) - t2im * x(iv,2,i,2,j)
          x2im =        t2re * x(iv,2,i,2,j) + t2im * x(iv,1,i,2,j)
          x3re = x2re - t3re * x(iv,1,i,3,j) + t3im * x(iv,2,i,3,j)
          x3im = x2im - t3re * x(iv,2,i,3,j) - t3im * x(iv,1,i,3,j)
          x0re = x1re - t4re * x(iv,1,i,4,j) + t4im * x(iv,2,i,4,j)
          x0im = x1im - t4re * x(iv,2,i,4,j) - t4im * x(iv,1,i,4,j)
          
          x1re =  2  * x1re -       x0re
          x1im =  2  * x1im -       x0im
          x4re =  2  * x2re -       x3re
          x4im =  2  * x2im -       x3im
          x2re =       x0re + C53 * x3re
          x2im =       x0im + C53 * x3im
          x3re = C53 * x0re -       x3re
          x3im = C53 * x0im -       x3im
          x0re =       x1re +       x4re
          x0im =       x1im +       x4im
          
          x1re = x1re          -       x4re
          x1im = x1im          -       x4im
          x4re = x(iv,1,i,0,j) - C51 * x0re
          x4im = x(iv,2,i,0,j) - C51 * x0im
          
          x1re =     x4re - C52 * x1re
          x1im =     x4im - C52 * x1im
          x4re = 2 * x4re -       x1re
          x4im = 2 * x4im -       x1im
          
          x(iv,1,i,0,j) =     x(iv,1,i,0,j) +       x0re
          x(iv,2,i,0,j) =     x(iv,2,i,0,j) +       x0im
          x(iv,1,i,3,j) =     x1re          - C54 * x3im
          x(iv,2,i,3,j) =     x1im          + C54 * x3re
          x(iv,1,i,2,j) = 2 * x1re          -       x(iv,1,i,3,j)
          x(iv,2,i,2,j) = 2 * x1im          -       x(iv,2,i,3,j)
          x(iv,1,i,4,j) =     x4re          - C54 * x2im
          x(iv,2,i,4,j) =     x4im          + C54 * x2re
          x(iv,1,i,1,j) = 2 * x4re          -       x(iv,1,i,4,j)
          x(iv,2,i,1,j) = 2 * x4im          -       x(iv,2,i,4,j)
        end do
      end do
    end do
    
  end procedure fxzm5a_c
  
end submodule fx5a
