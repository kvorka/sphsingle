submodule (fourier_transform) c2r
  implicit none ; contains
  
  module procedure fft_c2r_sub
    integer :: i
    
    call fxaddsub( m, x(1,1,0), x(1,2,0) )
    
    !$omp parallel do
    do i = 1, (this%n-2)/4
      call fxc2r( m, this%t(this%n+2*i-1), x(1,1,i), &
                                         & x(1,2,i), &
                                         & x(1,1,this%n/2-i), &
                                         & x(1,2,this%n/2-i)  )
    end do
    !$omp end parallel do
    
    if ( mod(this%n,4) == 0) then
      call fxrsc( m, +2, x(1,1,this%n/4) )
      call fxrsc( m, -2, x(1,2,this%n/4) )
    end if
    
    call fxztal( this%n, this%it(this%n/2-1), this%t, m, x )
    call fxzshf( this%n, this%it, 2*m, x )
    
  end procedure fft_c2r_sub
  
end submodule c2r