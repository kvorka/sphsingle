submodule (lateral_grid) init
  implicit none ; contains
  
  module procedure init_harmonics_sub
    
    call this%fft%init_sub( 3*jmax+1 )
    call this%lgp%init_sub( jmax, this%fft%n )
    
    this%nGrid = 2 * this%lgp%n * this%fft%n
    
  end procedure init_harmonics_sub
  
  module procedure deallocate_harmonics_sub
    
    call this%fft%deallocate_sub()
    call this%lgp%deallocate_sub()
    
  end procedure deallocate_harmonics_sub
  
end submodule init
