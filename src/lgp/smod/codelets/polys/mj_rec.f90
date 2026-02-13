submodule (lege_poly) poly_mj
  implicit none; contains
  
  module procedure mj_rec_sub
    integer                             :: i1, i2
    real(kind=dbl)                      :: cff1, cff2
    real(kind=dbl), pointer, contiguous :: cosx2(:,:)
    
    cff1 = this%fmj(1,ima)
    cff2 = this%fmj(2,ima)
    
    cosx2(1:ndbl,1:this%n_dbl) => this%cosx2
    
    do i2 = 1, this%n_dbl4, 4
      !$omp simd aligned (cosx2:alig)
      do i1 = 1, ndbl  
        pmj(i1,i2  ) = ( cff1 * cosx2(i1,i2  ) - cff2 ) * pmj1(i1,i2  ) - pmj(i1,i2  )
        pmj(i1,i2+1) = ( cff1 * cosx2(i1,i2+1) - cff2 ) * pmj1(i1,i2+1) - pmj(i1,i2+1)
        pmj(i1,i2+2) = ( cff1 * cosx2(i1,i2+2) - cff2 ) * pmj1(i1,i2+2) - pmj(i1,i2+2)
        pmj(i1,i2+3) = ( cff1 * cosx2(i1,i2+3) - cff2 ) * pmj1(i1,i2+3) - pmj(i1,i2+3)
      end do
    end do
    
    if ( this%n_dbl_div_4 ) then
      i2 = this%n_dbl4+1
      
      !$omp simd aligned (cosx2:alig)
      do i1 = 1, ndbl
        pmj(i1,i2  ) = ( cff1 * cosx2(i1,i2  ) - cff2 ) * pmj1(i1,i2  ) - pmj(i1,i2  )
        pmj(i1,i2+1) = ( cff1 * cosx2(i1,i2+1) - cff2 ) * pmj1(i1,i2+1) - pmj(i1,i2+1)
      end do
    end if
    
  end procedure mj_rec_sub
  
end submodule poly_mj