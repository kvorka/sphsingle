submodule (fourier_transform) fx4a
  implicit none; contains
  
  module procedure fxzm4a_c
    integer        :: i, j, iv
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im, t1re, t1im, t2re, t2im
    
    !$omp parallel private (i,iv,x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im, t1re, t1im, t2re, t2im)
    
    !$omp do
    do j = 0, k-1
      t1re = t(  6*j)
      t1im = t(1+6*j)
      t2re = t(2+6*j)
      t2im = t(3+6*j)
      
      do i = 1, l
        do iv = 1, howmany
          x0re = x(iv,1,i,0,j)
          x0im = x(iv,2,i,0,j)
          x1re = x(iv,1,i,1,j)
          x1im = x(iv,2,i,1,j)
          
          x2re = t2im * x(iv,2,i,2,j) - ( t2re * x(iv,1,i,2,j) - x0re )
          x2im = t2re * x(iv,2,i,2,j) + ( t2im * x(iv,1,i,2,j) - x0im )
          x3re = t2im * x(iv,2,i,3,j) - ( t2re * x(iv,1,i,3,j) - x1re )
          x3im = t2re * x(iv,2,i,3,j) + ( t2im * x(iv,1,i,3,j) - x1im )
          
          x0re = 2 * x0re - x2re
          x0im = 2 * x0im + x2im
          x1re = 2 * x1re - x3re
          x1im = 2 * x1im + x3im
          
          x(iv,1,i,2,j) = +( t1im * x1im - ( t1re * x1re - x0re ) )
          x(iv,2,i,2,j) = -( t1re * x1im + ( t1im * x1re - x0im ) )
          x(iv,1,i,0,j) = 2 * x0re - x(iv,1,i,2,j)
          x(iv,2,i,0,j) = 2 * x0im - x(iv,2,i,2,j)
          x(iv,1,i,1,j) = -( t1im * x3re - ( t1re * x3im + x2re ) )
          x(iv,2,i,1,j) = +( t1im * x3im + ( t1re * x3re - x2im ) )
          x(iv,1,i,3,j) =  2 * x2re - x(iv,1,i,1,j)
          x(iv,2,i,3,j) = -2 * x2im - x(iv,2,i,1,j)
        end do
      end do
    end do
    !$omp end parallel
    
  end procedure fxzm4a_c
  
end submodule fx4a
