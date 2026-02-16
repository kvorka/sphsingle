submodule (lege_poly) bwd_shuffle
  implicit none; contains
  
  module procedure bwd_shuffle_sub
    integer                             :: i1, i2
    real(kind=dbl)                      :: s11, s21, s31, s41, s12, s22, s32, s42, cx1, cx2
    real(kind=dbl), pointer, contiguous :: cosx(:,:,:)
    
    cosx(1:ndbl,1:2,1:this%n_dbl_2) => this%cosx
    
    do i2 = 1, this%n_dbl_2
      !$omp simd aligned (cosx:alig)
      do i1 = 1, ndbl
        cx1 = cosx(i1,1,i2)
        cx2 = cosx(i1,2,i2)
        
        s11 = swork(i1,1,i2)
        s21 = swork(i1,3,i2)
        s31 = swork(i1,2,i2)
        s41 = swork(i1,4,i2)
        s12 = swork(i1,5,i2)
        s22 = swork(i1,7,i2)
        s32 = swork(i1,6,i2)
        s42 = swork(i1,8,i2)
        
        grid(i1,1,i2,1) = s21 * cx1 + s11
        grid(i1,1,i2,2) = s21 * cx1 - s11
        grid(i1,1,i2,3) = s41 * cx1 + s31
        grid(i1,1,i2,4) = s41 * cx1 - s31
        grid(i1,2,i2,1) = s22 * cx2 + s12
        grid(i1,2,i2,2) = s22 * cx2 - s12
        grid(i1,2,i2,3) = s42 * cx2 + s32
        grid(i1,2,i2,4) = s42 * cx2 - s32
      end do
    end do
    
  end procedure bwd_shuffle_sub

end submodule bwd_shuffle