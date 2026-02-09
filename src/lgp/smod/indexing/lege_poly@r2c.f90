submodule (lege_poly) r2c
  implicit none; contains
  
  module procedure index_fwd_sub
    integer                     :: m, j, jm, mj
    real(kind=dbl), allocatable :: cab(:,:)
    
    allocate( cab(2,this%jms) )
    
    call this%is_rescale_sub( rcab )
    call this%fwd_indx_sub( rcab, cab )
    
    do j = 0, this%jmax
      !m = 0
        jm = j*(j+1)/2+1
        cjm(jm) = cmplx( cab(1,j+1), 0._dbl, kind=dbl )
      
      !$omp simd
      do m = 1, j
        mj = m*(this%jmax+1)-m*(m+1)/2+j+1
        cjm(jm+m) = cmplx( cab(1,mj), cab(2,mj), kind=dbl )
      end do
    end do
    
    deallocate( cab )
    
  end procedure index_fwd_sub
  
end submodule r2c
