submodule (lege_poly) poly_mj
  implicit none; contains
  
  module procedure mj_rec_sub
    integer                             :: i1, i2
    real(kind=dbl)                      :: cff1, cff2
    real(kind=dbl), pointer, contiguous :: cosx2(:,:)
    
    cosx2(1:16,1:this%n_16) => this%cosx2
    
    cff1 = this%fmj(1,ima)
    cff2 = this%fmj(2,ima)
    
    do i2 = 1, this%n_16
      !$omp simd aligned (cosx2:alig)
      do i1 = 1, 16
        pmj(i1,i2) = ( cff1 * cosx2(i1,i2) - cff2 ) * pmj1(i1,i2) - pmj(i1,i2)
      end do
    end do
    
  end procedure mj_rec_sub
  
end submodule poly_mj