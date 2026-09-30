submodule (lege_poly) r2c
  implicit none; contains
  
  module procedure index_fwd_sub
    real(kind=dbl), allocatable :: cab(:)
    
    !! First rescale required by the Ishioka recursion. This call is not
    !! inlined because the same rescale is required before the backword
    !! transform. Vectorized inside.
    call is_rescale_sub( this%nrma, this%amj, rcab )
    
    !! Reindexing after transform. Quadruplets of rescaled coefficients
    !! are synthethysed. These include real/imaginary, odd degree/even degree
    !! components, respectively. Afterwards, reindexing from order-fast mj to 
    !! degree-fast jm indexing and synthethysing the real/imaginary parts into cmplx.
    allocate( cab(2*this%jms) )
    
    call fwd_rxd_sub( this%jmax, this%emj, rcab, cab )
    call fwd_r2c_sub( this%jmax, cab, cjm )
    
    deallocate( cab )
    
  end procedure index_fwd_sub
  
end submodule r2c
