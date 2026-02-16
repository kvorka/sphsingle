submodule (lege_poly) fwd_sum2
  implicit none; contains
  
  module procedure fwd_sum2_sub
    integer                             :: i1, i2
    real(kind=dbl)                      :: cff1, cff2
    real(kind=dbl), pointer, contiguous :: cosx2(:,:,:)
    
    cff1 = this%fmj(1,ima)
    cff2 = this%fmj(2,ima)
    
    cosx2(1:ndbl,1:2,1:this%n_dbl_2) => this%cosx2
    
    !$omp simd
    do i1 = 1, ndbl
      acc(i1,1) = zero
      acc(i1,2) = zero
      acc(i1,3) = zero
      acc(i1,4) = zero
      
      acc2(i1,1) = zero
      acc2(i1,2) = zero
      acc2(i1,3) = zero
      acc2(i1,4) = zero
    end do
    
    do i2 = 1, this%n_dbl_2
      !$omp simd aligned (cosx2:alig)
      do i1 = 1, ndbl
        pmj(i1,1,i2) = ( cff1 * cosx2(i1,1,i2) - cff2 ) * pmj1(i1,1,i2) - pmj(i1,1,i2)
        pmj(i1,2,i2) = ( cff1 * cosx2(i1,2,i2) - cff2 ) * pmj1(i1,2,i2) - pmj(i1,2,i2)
        
        acc(i1,1) = acc(i1,1) + pmj(i1,1,i2) * swork(i1,1,i2)
        acc(i1,2) = acc(i1,2) + pmj(i1,1,i2) * swork(i1,2,i2)
        acc(i1,3) = acc(i1,3) + pmj(i1,1,i2) * swork(i1,3,i2)
        acc(i1,4) = acc(i1,4) + pmj(i1,1,i2) * swork(i1,4,i2)
        
        acc2(i1,1) = acc2(i1,1) + pmj(i1,2,i2) * swork(i1,5,i2)
        acc2(i1,2) = acc2(i1,2) + pmj(i1,2,i2) * swork(i1,6,i2)
        acc2(i1,3) = acc2(i1,3) + pmj(i1,2,i2) * swork(i1,7,i2)
        acc2(i1,4) = acc2(i1,4) + pmj(i1,2,i2) * swork(i1,8,i2)
      end do
    end do
    
    !$omp simd
    do i1 = 1, ndbl
      acc(i1,1) = acc(i1,1) + acc2(i1,1)
      acc(i1,2) = acc(i1,2) + acc2(i1,2)
      acc(i1,3) = acc(i1,3) + acc2(i1,3)
      acc(i1,4) = acc(i1,4) + acc2(i1,4)
      
      cr(1) = cr(1) + acc(i1,1)
      cr(2) = cr(2) + acc(i1,2)
      cr(3) = cr(3) + acc(i1,3)
      cr(4) = cr(4) + acc(i1,4)
    end do
    
  end procedure fwd_sum2_sub
  
end submodule fwd_sum2