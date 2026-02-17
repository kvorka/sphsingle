submodule (lateral_grid) init
  implicit none ; contains
  
  module procedure init_harmonics_sub
    integer :: nL, nF
    
    nL = (3*jmax/2+1)/2+ndbl+1-mod((3*jmax/2+1)/2+1,ndbl)
    nF = prime_adjustement_sub(3*jmax+1)
    
    this%nGrid = 2 * nL * nF
    
    call this%fft%init_sub( nF )
    call this%lgp%init_sub( jmax, nL, real(nF, kind=dbl) )
    
  end procedure init_harmonics_sub
  
  module procedure deallocate_harmonics_sub
    
    call this%fft%deallocate_sub()
    call this%lgp%deallocate_sub()
    
  end procedure deallocate_harmonics_sub
  
end submodule init
