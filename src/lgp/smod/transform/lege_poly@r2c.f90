submodule (lege_poly) r2c
  implicit none; contains
  
  module procedure index_fwd_sub
    real(kind=dbl), allocatable :: cab(:)
    
    !! Reindexing after transform. Quadruplets of coefficients are rescaled and synthethysed. 
    !! These include real/imaginary, odd degree/even degree components, respectively. Afterwards, 
    !! reindexing from order-fast mj to degree-fast jm indexing and casting the real/imaginary 
    !! parts into cmplx is carried out.
    allocate( cab(2*this%jms) )
    
    call fwd_rxd_sub( this%jmax, this%emj, this%amj, rcab, cab )
    call fwd_r2c_sub( this%jmax, cab, cjm )
    
    deallocate( cab )
    
  end procedure index_fwd_sub
  
end submodule r2c
