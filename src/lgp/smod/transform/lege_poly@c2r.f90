submodule (lege_poly) c2r
  implicit none; contains
  
  module procedure index_bwd_sub
    real(kind=dbl), allocatable :: cab(:)
    
    !! Reindexing of real/imaginary parts of the original sequence from order-fast jm 
    !! to degree-fast mj indexing. Afterwards, reindexing and rescaling for transform.
    !! Quadruplets of the rescaled coefficietns order real/imaginary, odd/even degrees
    !! components, respectively, for cache friendly behaviour.
    allocate( cab(2*this%jms) )
    
    call bwd_c2r_sub( this%jmax, cjm, cab )
    call bwd_rxd_sub( this%jmax, this%emj, this%amj, cab, rcab )
    
    deallocate( cab )
    
  end procedure index_bwd_sub
  
end submodule c2r