submodule (fourier_transform) c2r
  implicit none ; contains
  
  module procedure fft_c2r_sub
    integer        :: i, iv
    real(kind=dbl) :: addre, addim, subre, subim, t1, t2
    
    do iv = 1, howmany
      addre     =             x(iv,1,0)
      x(iv,1,0) = x(iv,1,0) + x(iv,2,0)
      x(iv,2,0) = addre     - x(iv,2,0)
    end do
    
    do i = 1, (this%n-2)/4
      t1 = this%t(this%n+2*i-1)
      t2 = this%t(this%n+2*i  )
      
      do iv = 1, howmany
        addre = x(iv,1,i) + x(iv,1,this%n/2-i)
        subre = x(iv,1,i) - x(iv,1,this%n/2-i)
        addim = x(iv,2,i) + x(iv,2,this%n/2-i)
        subim = x(iv,2,i) - x(iv,2,this%n/2-i)
        
        x(iv,1,i) = addre - subre * t2 - addim * t1
        x(iv,2,i) = subim - addim * t2 + subre * t1
        
        x(iv,1,this%n/2-i) = -x(iv,1,i) + 2 * addre
        x(iv,2,this%n/2-i) = +x(iv,2,i) - 2 * subim
      end do
    end do
    
    if ( mod(this%n,4) == 0) then
      do iv = 1, howmany
        x(iv,1,this%n/4) = +x(iv,1,this%n/4) * 2
        x(iv,2,this%n/4) = -x(iv,2,this%n/4) * 2
      end do
    end if
    
    call this%fxztal( howmany, x )
    call this%fxzshf( howmany, x )
    
  end procedure fft_c2r_sub
  
end submodule c2r