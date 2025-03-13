submodule (lege_poly) c2r
  implicit none; contains
  
  module procedure index_bwd_sub
    integer                             :: m, j, jm, mj
    type(c_ptr)                         :: c_cab
    real(kind=dbl), pointer, contiguous :: cab(:,:)
    
    call alloc_aligned_2d_sub( 2, this%jms, c_cab, cab )
    
    do m = 0, this%jmax
      do j = m, this%jmax
        mj = m*(this%jmax+1)-m*(m+1)/2+j+1
        jm = j*(j+1)/2+m+1
        
        cab(1,mj) = cjm(jm)%re
        cab(2,mj) = cjm(jm)%im
      end do
    end do
    
    call bwd_indx_c( this%jmax, this%emj, cab, rcab )
    call is_rescale_c( this%nrma, this%amj, rcab )
    
    call free( c_cab )
    
  end procedure index_bwd_sub
  
end submodule c2r