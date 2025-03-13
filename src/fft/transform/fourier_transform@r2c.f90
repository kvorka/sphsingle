submodule (fourier_transform) r2c
  implicit none; contains
  
  module procedure fft_r2c_sub
    integer        :: i, iv, iv1
    real(kind=dbl) :: addre, addim, subre, subim, t1, t2
    
    call this%fxztal( howmany, x )
    call this%fxzshf( howmany, x )
    
    do iv = 1, howmany, 16
      !$omp simd
      do iv1 = 0, 15
        addre         =                 x(iv1+iv,1,0)
        x(iv1+iv,1,0) = x(iv1+iv,1,0) + x(iv1+iv,2,0)
        x(iv1+iv,2,0) = addre         - x(iv1+iv,2,0)
      end do
    end do
    
    do i = 1, (this%n-2)/4
      t1 = this%t(this%n+2*i-1)
      t2 = this%t(this%n+2*i  )
      
      do iv = 1, howmany, 16
        !$omp simd
        do iv1 = 0, 15
          addre = -x(iv1+iv,1,i) - x(iv1+iv,1,this%n/2-i)
          subre = +x(iv1+iv,1,i) - x(iv1+iv,1,this%n/2-i)
          addim = +x(iv1+iv,2,i) + x(iv1+iv,2,this%n/2-i)
          subim = -x(iv1+iv,2,i) + x(iv1+iv,2,this%n/2-i)
          
          x(iv1+iv,1,i) = -( addre - subre * t2 - addim * t1 ) / 2
          x(iv1+iv,2,i) = +( subim - addim * t2 + subre * t1 ) / 2
          
          x(iv1+iv,1,this%n/2-i) = -x(iv1+iv,1,i) - addre
          x(iv1+iv,2,this%n/2-i) = +x(iv1+iv,2,i) - subim
        end do
      end do
    end do
    
    if ( mod(this%n,4) == 0) then
      do iv = 1, howmany, 16
        !$omp simd
        do iv1 = 0, 15
          x(iv1+iv,1,this%n/4) = +x(iv1+iv,1,this%n/4)
          x(iv1+iv,2,this%n/4) = -x(iv1+iv,2,this%n/4)
        end do
      end do
    end if
    
  end procedure fft_r2c_sub
  
end submodule r2c