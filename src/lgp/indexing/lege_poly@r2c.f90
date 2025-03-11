submodule (lege_poly) r2c
  implicit none; contains
  
  module procedure index_fwd_sub
    integer                 :: m, j, jm, mj
    type(c_ptr)             :: c_arr
    real(kind=dbl), pointer :: cab(:,:)
    
    c_arr = malloc( 32, 2 * this%jms * size_d )
    call c_f_pointer( c_arr, cab, [2,this%jms] )
    
    call is_rescale_c( this%nrma, this%amj, rcab )
    call fwd_indx_c( this%jmax, this%emj, rcab, cab )
    
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
    
    call free( c_arr )
    
  end procedure index_fwd_sub
  
end submodule r2c
