submodule (fourier_transform) fx3a
  implicit none; contains
  
  module procedure fxzm3a
    integer        :: i1, i2, i3, i4
    real(kind=dbl) :: t1re, t1im, t2re, t2im, x0re, x0im, x1re, x1im, x2re, x2im, x01, x02
    
    !$omp parallel do private (i1,i2,i3,t1re,t1im,t2re,t2im,x0re,x0im,x1re,x1im,x2re,x2im,x01,x02)
    do i4 = 0, k-1
      i1 = 2 * i4
      
      t1re = t(1,i1  )
      t1im = t(2,i1  )            
      t2re = t(1,i1+1)
      t2im = t(2,i1+1)
      
      do i3 = 1, l/3
        do i2 = 1, m
          !$omp simd
          do i1 = 1, ndbl
            x0re = x(i1,i2,1,i3,1,i4)
            x0im = x(i1,i2,2,i3,1,i4)
            x2re = x(i1,i2,1,i3,2,i4)
            x2im = x(i1,i2,2,i3,2,i4)
            
            x01  = t1im * x0im
            x02  = t1im * x0re
            x1re = t2re * x2re
            x1im = t2re * x2im
            
            x01  = +t1re * x0re - x01
            x02  = +t1re * x0im + x02
            x1re = -t2im * x2im + x1re
            x1im = +t2im * x2re + x1im
            
            x1re = x01 - x1re
            x1im = x02 - x1im
            x2re = x01 + x01
            x2im = x02 + x02
            
            x01 = x2re - x1re
            x02 = x2im - x1im
            
            x0re = x(i1,i2,1,i3,0,i4)
            x0im = x(i1,i2,2,i3,0,i4)
            
            x2re = x0re + C31 * x01
            x2im = x0im + C31 * x02
            
            x0re = x01 + x0re
            x0im = x02 + x0im
            
            x(i1,i2,1,i3,0,i4) = x0re
            x(i1,i2,2,i3,0,i4) = x0im
            
            x1im = x2re + C32 * x1im
            x1re = x2im - C32 * x1re
            
            x2re = x2re + x2re
            x2im = x2im + x2im
            
            x(i1,i2,1,i3,2,i4) = x1im
            x(i1,i2,2,i3,2,i4) = x1re
            
            x2re = x2re - x1im
            x2im = x2im - x1re
            
            x(i1,i2,1,i3,1,i4) = x2re
            x(i1,i2,2,i3,1,i4) = x2im
          end do
        end do
      end do
    end do
    !$omp end parallel do
    
  end procedure fxzm3a
  
end submodule fx3a
