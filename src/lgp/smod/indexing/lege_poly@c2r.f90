submodule (lege_poly) c2r
  implicit none; contains
  
  module procedure index_bwd_sub
    integer                     :: m, j, jm, mj
    real(kind=dbl), allocatable :: cab(:,:)
    
    allocate( cab(2,this%jms) )
    
    do m = 0, this%jmax
      mj = m*(this%jmax+1)-m*(m+1)/2+1
      
      !$omp simd
      do j = m, this%jmax
        jm = j*(j+1)/2+m+1
        
        cab(1,mj+j) = cjm(jm)%re
        cab(2,mj+j) = cjm(jm)%im
      end do
    end do
    
    call bwd_indx_sub( this%jmax, this%emj, cab, rcab )
    call is_rescale_sub( this%nrma, this%amj, rcab )
    
    deallocate( cab )
    
  end procedure index_bwd_sub
  
end submodule c2r