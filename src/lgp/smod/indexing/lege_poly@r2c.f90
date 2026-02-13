submodule (lege_poly) r2c
  implicit none; contains
  
  module procedure index_fwd_sub
    integer                     :: m, j, jm, mj
    real(kind=dbl), allocatable :: cab(:,:)
    
    allocate( cab(2,this%jms) )
    
    call this%is_rescale_sub( rcab )
    call this%fwd_indx_sub( rcab, cab )
    
    m = 0
      !$omp simd
      do j = 0, this%jmax
        jm = j*(j+1)/2+1
        
        cjm(jm)%re = cab(1,1+j)
        cjm(jm)%im = zero
      end do
    
    do m = 1, this%jmax
      mj = m*(this%jmax+1)-m*(m+1)/2+1
      
      !$omp simd
      do j = m, this%jmax
        jm = j*(j+1)/2+m+1
        
        cjm(jm)%re = cab(1,mj+j)
        cjm(jm)%im = cab(2,mj+j)
      end do
    end do
    
    deallocate( cab )
    
  end procedure index_fwd_sub
  
end submodule r2c
