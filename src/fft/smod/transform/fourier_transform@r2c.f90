submodule (fourier_transform) r2c
  implicit none; contains
  
  module procedure fft_r2c_sub
    integer :: i
    
    call fxztal( this%n, this%it, this%t, m, x )
    
    call fxrc0( m, x(1,1,0), x(1,2,0) )
    
    do i = 1, (this%n-2)/4
      call fxr2c( m, this%t(this%n+2*i-1), x(1,1,i), &
                                         & x(1,2,i), &
                                         & x(1,1,this%n/2-i), &
                                         & x(1,2,this%n/2-i)  )
    end do
    
    if ( mod(this%n,4) == 0 ) call fxrsc( m, -1._dbl, x(1,2,this%n/4) )
    
  end procedure fft_r2c_sub
  
end submodule r2c