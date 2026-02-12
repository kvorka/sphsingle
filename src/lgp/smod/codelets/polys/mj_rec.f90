submodule (lege_poly) poly_mj
  implicit none; contains
  
  module procedure mj_rec_sub
    integer                             :: i1, i2
    real(kind=dbl)                      :: cff1, cff2
    real(kind=dbl), pointer, contiguous :: cosx2(:,:,:)
    
    cff1 = this%fmj(1,ima)
    cff2 = this%fmj(2,ima)
    
    cosx2(1:ndbl,1:2,1:this%n_step) => this%cosx2
    
    do i2 = 1, this%nUnroll3, 3
      !$omp simd aligned (cosx2:alig)
      do i1 = 1, ndbl
        pmj(i1,1,i2  ) = ( cff1 * cosx2(i1,1,i2  ) - cff2 ) * pmj1(i1,1,i2  ) - pmj(i1,1,i2  )
        pmj(i1,2,i2  ) = ( cff1 * cosx2(i1,2,i2  ) - cff2 ) * pmj1(i1,2,i2  ) - pmj(i1,2,i2  )
        pmj(i1,1,i2+1) = ( cff1 * cosx2(i1,1,i2+1) - cff2 ) * pmj1(i1,1,i2+1) - pmj(i1,1,i2+1)
        pmj(i1,2,i2+1) = ( cff1 * cosx2(i1,2,i2+1) - cff2 ) * pmj1(i1,2,i2+1) - pmj(i1,2,i2+1)
        pmj(i1,1,i2+2) = ( cff1 * cosx2(i1,1,i2+2) - cff2 ) * pmj1(i1,1,i2+2) - pmj(i1,1,i2+2)
        pmj(i1,2,i2+2) = ( cff1 * cosx2(i1,2,i2+2) - cff2 ) * pmj1(i1,2,i2+2) - pmj(i1,2,i2+2)
      end do
    end do
    
    do i2 = this%nUnroll3+1, this%n_step
      !$omp simd aligned (cosx2:alig)
      do i1 = 1, ndbl
        pmj(i1,1,i2) = ( cff1 * cosx2(i1,1,i2) - cff2 ) * pmj1(i1,1,i2) - pmj(i1,1,i2)
        pmj(i1,2,i2) = ( cff1 * cosx2(i1,2,i2) - cff2 ) * pmj1(i1,2,i2) - pmj(i1,2,i2)
      end do
    end do
    
  end procedure mj_rec_sub
  
end submodule poly_mj