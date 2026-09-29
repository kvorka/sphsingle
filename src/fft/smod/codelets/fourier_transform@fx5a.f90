submodule (fourier_transform) fx5a
  implicit none; contains
  
  module procedure fxzm5a
    integer        :: i1, i2, i3, i4
    real(kind=dbl) :: t1re, t1im, t2re, t2im, t3re, t3im, t4re, t4im, &
                      x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im, &
                      x4re, x4im, x01,  x02,  x03,  x04
    
    !$omp parallel do private (i1,i2,i3,t1re,t1im,t2re,t2im,t3re,t3im,t4re,t4im, &
    !$omp& x0re,x0im,x1re,x1im,x2re,x2im,x3re,x3im,x4re,x4im,x01,x02,x03,x04)
    do i4 = 0, k-1
      i1 = 4 * i4
      
      t1re = t(1,i1  )
      t1im = t(2,i1  )
      t2re = t(1,i1+1)
      t2im = t(2,i1+1)
      t3re = t(1,i1+2)
      t3im = t(2,i1+2)
      t4re = t(1,i1+3)
      t4im = t(2,i1+3)
      
      do i3 = 1, l/5
        do i2 = 1, m
          !$omp simd
          do i1 = 1, ndbl
            x0re = x(i1,i2,1,i3,1,i4)
            x0im = x(i1,i2,2,i3,1,i4)
            
            x1re = t1re * x0re
            x1im = t1re * x0im
            
            x3re = x(i1,i2,1,i3,2,i4)
            x3im = x(i1,i2,2,i3,2,i4)
            
            x2re = t2re * x3re
            x2im = t2re * x3im
            
            x1re = x1re - t1im * x0im
            x1im = x1im + t1im * x0re
            x2re = x2re - t2im * x3im
            x2im = x2im + t2im * x3re
            
            x01 = x(i1,i2,1,i3,3,i4)
            x02 = x(i1,i2,2,i3,3,i4)
            
            x3re = x2re - t3re * x01
            x3im = x2im - t3re * x02
            
            x03 = x(i1,i2,1,i3,4,i4)
            x04 = x(i1,i2,2,i3,4,i4)
            
            x0re = x1re - t4re * x03
            x0im = x1im - t4re * x04
            
            x3re = x3re + t3im * x02
            x3im = x3im - t3im * x01
            x0re = x0re + t4im * x04
            x0im = x0im - t4im * x03
            
            x1re = 2 * x1re - x0re
            x1im = 2 * x1im - x0im
            x4re = 2 * x2re - x3re
            x4im = 2 * x2im - x3im
            
            x2re = x0re + C53 * x3re
            x2im = x0im + C53 * x3im
            
            x3re = C53 * x0re - x3re
            x3im = C53 * x0im - x3im
            
            x0re = x1re + x4re
            x0im = x1im + x4im
            
            x1re = x1re - x4re
            x1im = x1im - x4im
            
            x4re = x(i1,i2,1,i3,0,i4)
            x4im = x(i1,i2,2,i3,0,i4)
            
            x01 = x4re + x0re
            x02 = x4im + x0im
            
            x(i1,i2,1,i3,0,i4) = x01
            x(i1,i2,2,i3,0,i4) = x02
            
            x4re = x4re - C51 * x0re
            x4im = x4im - C51 * x0im
            
            x1re = x4re - C52 * x1re
            x1im = x4im - C52 * x1im
            
            x4re = 2 * x4re - x1re
            x4im = 2 * x4im - x1im
            x3im = x1re - C54 * x3im
            x3re = x1im + C54 * x3re
            
            x(i1,i2,1,i3,3,i4) = x3im
            x(i1,i2,2,i3,3,i4) = x3re
            
            x1re = 2 * x1re - x3im
            x1im = 2 * x1im - x3re
            x2im = x4re - C54 * x2im
            x2re = x4im + C54 * x2re
            
            x(i1,i2,1,i3,2,i4) = x1re
            x(i1,i2,2,i3,2,i4) = x1im
            
            x4re = 2 * x4re - x2im
            x4im = 2 * x4im - x2re
            
            x(i1,i2,1,i3,4,i4) = x2im
            x(i1,i2,2,i3,4,i4) = x2re
            x(i1,i2,1,i3,1,i4) = x4re
            x(i1,i2,2,i3,1,i4) = x4im
          end do
        end do
      end do
    end do
    !$omp end parallel do
    
  end procedure fxzm5a
  
end submodule fx5a
