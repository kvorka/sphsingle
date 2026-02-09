submodule (lege_poly) r2c
  implicit none; contains
  
  module procedure index_fwd_sub
    integer                             :: m, j, jm, mj
    type(c_ptr)                         :: c_cab
    real(kind=dbl), pointer, contiguous :: cab(:,:)
    
    call alloc_aligned_2d_sub( 2, this%jms, c_cab, cab )
    
    call is_rescale_c( this%nrma, this%amj, rcab )
    call this%fwd_indx_sub( rcab, cab )
    
    do j = 0, this%jmax
      m = 0
        mj = m*(this%jmax+1)-m*(m+1)/2+j+1
        jm = j*(j+1)/2+m+1
        
        cjm(jm) = cmplx( cab(1,mj), 0._dbl, kind=dbl )
      
      do m = 1, j
        mj = m*(this%jmax+1)-m*(m+1)/2+j+1
        jm = j*(j+1)/2+m+1
        
        cjm(jm) = cmplx( cab(1,mj), cab(2,mj), kind=dbl )
      end do
    end do
    
    call free( c_cab )
    
  end procedure index_fwd_sub
  
end submodule r2c
