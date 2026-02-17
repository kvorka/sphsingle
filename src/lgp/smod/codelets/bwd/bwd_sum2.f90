submodule (lege_poly) bwd_sum2
  implicit none; contains
  
  module procedure bwd_sum2_sub
    integer                             :: i1, i2
    real(kind=dbl)                      :: p1, p2, c1, c2, c3, c4, cff1, cff2
    real(kind=dbl), pointer, contiguous :: cosx2(:,:,:)
    
    c1 = cc(1)
    c2 = cc(2)
    c3 = cc(3)
    c4 = cc(4)
    
    cff1 = this%fmj(1,ma)
    cff2 = this%fmj(2,ma)
    
    cosx2(1:ndbl,1:2,1:this%n_dbl_2) => this%cosx2
    
    do i2 = 1, this%n_dbl_2
      !$omp simd aligned (cosx2:alig)
      do i1 = 1, ndbl
        pmj(i1,1,i2) = ( cff1 * cosx2(i1,1,i2) - cff2 ) * pmj1(i1,1,i2) - pmj(i1,1,i2)
        pmj(i1,2,i2) = ( cff1 * cosx2(i1,2,i2) - cff2 ) * pmj1(i1,2,i2) - pmj(i1,2,i2)
        
        swork(i1,1,i2) = swork(i1,1,i2) + pmj(i1,1,i2) * c1
        swork(i1,2,i2) = swork(i1,2,i2) + pmj(i1,1,i2) * c2
        swork(i1,3,i2) = swork(i1,3,i2) + pmj(i1,1,i2) * c3
        swork(i1,4,i2) = swork(i1,4,i2) + pmj(i1,1,i2) * c4
        swork(i1,5,i2) = swork(i1,5,i2) + pmj(i1,2,i2) * c1
        swork(i1,6,i2) = swork(i1,6,i2) + pmj(i1,2,i2) * c2
        swork(i1,7,i2) = swork(i1,7,i2) + pmj(i1,2,i2) * c3
        swork(i1,8,i2) = swork(i1,8,i2) + pmj(i1,2,i2) * c4
      end do
    end do
    
  end procedure bwd_sum2_sub
  
end submodule bwd_sum2