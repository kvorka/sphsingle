submodule (lege_poly) init
  implicit none; contains
  
  module procedure init_lege_sub
    
    this%jmax = jmax
    this%jms  = jmax*(jmax+1)/2+jmax+1
    
    this%n     = n
    this%n_dbl = n / ndbl
    
    this%n_dbl4 = ( this%n_dbl / 4 ) * 4
    this%n_dbl_div_4 = .not. ( this%n_dbl == this%n_dbl4 )
    
    this%nFreq = 4 * this%n * ( this%jmax+1 )
    
    call this%get_nma_sub()
    call this%roots_sub()
    call this%coeffs_sub()
    call this%pmm_sub()
    
    this%wght = this%wght / wfac
    
  end procedure init_lege_sub
  
  module procedure deallocate_lege_sub
    
    call free( this%c_cosx  )
    call free( this%c_pmm   )
    call free( this%c_cosx2 )
    call free( this%c_wght  )
    
    this%cosx  => null()
    this%pmm   => null()
    this%cosx2 => null()
    this%wght  => null()
    
    deallocate( this%amj )
    deallocate( this%emj )
    deallocate( this%fmj )
    deallocate( this%mamj )
    
  end procedure deallocate_lege_sub
  
end submodule init