submodule (lege_poly) init
  implicit none; contains
  
  module procedure init_lege_sub
    integer :: j, m
    
    this%n = n
    this%jmax  = jmax
    this%jms   = jmax*(jmax+1)/2+jmax+1
    
    this%nrma = 0
      do m = 0, this%jmax
        this%nrma = this%nrma+1
        
        if ( m < this%jmax) then
          do j = 1, (this%jmax-1-m)/2
            this%nrma = this%nrma+1
          end do
          
          this%nrma = this%nrma+1
        end if
      end do
    
    call this%roots_sub()
    call this%coeffs_sub()
    
    this%rw(:,4) = this%rw(:,4) / wfac
    
  end procedure init_lege_sub
  
  module procedure deallocate_lege_sub
    
    if ( allocated(this%rw)  ) deallocate( this%rw  ) 
    if ( allocated(this%emj) ) deallocate( this%emj )
    if ( allocated(this%fmj) ) deallocate( this%fmj )
    
  end procedure deallocate_lege_sub
  
end submodule init