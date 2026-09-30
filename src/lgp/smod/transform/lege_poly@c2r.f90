submodule (lege_poly) c2r
  implicit none; contains
  
  module procedure index_bwd_sub
    real(kind=dbl), allocatable :: cab(:)
    
    !! Reindexing of real/imaginary parts of the original sequence 
    !! from order-fast jm to degree-fast mj indexing. Afterwards,
    !! reindexing for transform. Quadruplets of rescaled coefficients
    !! are prepared. These include real/imaginary, odd degree/even degree
    !! components, respectively, for easier caching.
    allocate( cab(2*this%jms) )
    
    call bwd_c2r_sub( this%jmax, cjm, cab )
    call bwd_rxd_sub( this%jmax, this%emj, cab, rcab )
    
    deallocate( cab )
    
    !! Last rescale required by the Ishioka recursion. This call is not
    !! inlined because the same rescale is required after the forward
    !! transform. Vectorized inside.
    call is_rescale_sub( this%nrma, this%amj, rcab )
    
  end procedure index_bwd_sub
  
end submodule c2r