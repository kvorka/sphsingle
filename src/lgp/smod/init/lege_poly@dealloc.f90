submodule (lege_poly) dealloc
  implicit none; contains
  
  module procedure deallocate_lege_sub
    
    call free( this%c_cosx  )
    call free( this%c_sinx  )
    call free( this%c_cosx2 )
    call free( this%c_wght  )
    
    this%cosx  => null()
    this%sinx  => null()
    this%cosx2 => null()
    this%wght  => null()
    
    deallocate( this%amj  )
    deallocate( this%emj  )
    deallocate( this%fmj  )
    deallocate( this%mamj )
    
  end procedure deallocate_lege_sub
  
end submodule dealloc