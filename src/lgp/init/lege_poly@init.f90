submodule (lege_poly) init
  implicit none; contains
  
  module procedure init_lege_sub
    integer :: j, m
    
    this%n = n
    this%jmax  = jmax
    this%jms   = jmax*(jmax+1)/2+jmax+1
    
    allocate( this%mma(0:this%jmax) )
    
    this%nrma = 0
      do m = 0, this%jmax
        this%nrma   = this%nrma+1
        this%mma(m) = this%nrma
        
        do j = 1, (this%jmax-m)/2
          this%nrma = this%nrma+1
        end do
        
        if ( mod((this%jmax-m),2) /= 0 ) then
          this%nrma = this%nrma+1
        end if
      end do
    
    call this%roots_sub()
    call this%coeffs_sub()
    call this%pmm_sub()
    
    this%wght = this%wght / wfac
    
  end procedure init_lege_sub
  
  module procedure deallocate_lege_sub
    
    if ( c_associated(this%c_cosx)  ) call free( this%c_cosx  )
    if ( c_associated(this%c_pmm)   ) call free( this%c_pmm   )
    if ( c_associated(this%c_cosx2) ) call free( this%c_cosx2 )
    if ( c_associated(this%c_wght)  ) call free( this%c_wght  )
    
    this%cosx  => null()
    this%pmm   => null()
    this%cosx2 => null()
    this%wght  => null()
    
    if ( allocated(this%amj) ) deallocate( this%emj )
    if ( allocated(this%emj) ) deallocate( this%emj )
    if ( allocated(this%fmj) ) deallocate( this%fmj )
    if ( allocated(this%mma) ) deallocate( this%mma )
    
  end procedure deallocate_lege_sub
  
end submodule init