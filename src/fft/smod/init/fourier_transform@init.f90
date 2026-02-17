submodule (fourier_transform) init
  implicit none; contains
  
  module procedure fft_init_sub
    integer :: i
    
    this%n = prime_adjustement_sub(n)
    
    allocate( this%it(this%n/2)  ) ; this%it = 0
    allocate( this%t(3*this%n/2) ) ; this%t = zero
    
    call fxzini( this%n/2, this%it, this%t )
    
    do i = 1, (this%n-2)/4
      this%t(this%n+2*i-1) = cos(2 * pi * i / this%n)
      this%t(this%n+2*i  ) = sin(2 * pi * i / this%n)        
    end do
    
  end procedure fft_init_sub
  
  module procedure fft_deallocate_sub
    
    deallocate( this%it )
    deallocate( this%t  )
    
  end procedure fft_deallocate_sub
  
end submodule init