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
    
    call this%bwd_indx_sub( cab, rcab )
    call this%is_rescale_sub( rcab )
    
    deallocate( cab )
    
  end procedure index_bwd_sub
  
end submodule c2r