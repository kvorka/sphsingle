submodule (fourier_transform) dealloc
  implicit none; contains
  
  module procedure fft_deallocate_sub
    
    deallocate( this%it )
    deallocate( this%t  )
    
  end procedure fft_deallocate_sub
  
end submodule dealloc