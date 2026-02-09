submodule (lege_poly) poly_mj
  implicit none; contains
  
  module procedure mj_rec_sub
    integer                             :: i1, i2, n64
    real(kind=dbl)                      :: cff0, cff1, cff2
    real(kind=dbl), pointer, contiguous :: p_pmj(:), p_pmj1(:), p_cosx2(:)
    
    n64 = (this%n/64)*64
    
    cff1 = this%fmj(1,ima)
    cff2 = this%fmj(2,ima)
    
    do i2 = 1, n64, 64
      p_cosx2(1:64) => this%cosx2(i2:i2+63)
      p_pmj(1:64)   => pmj(i2:i2+63)
      p_pmj1(1:64)  => pmj1(i2:i2+63)
      
      !$omp simd aligned (p_pmj,p_pmj1,p_cosx2:alig)
      do i1 = 1, 64
        p_pmj(i1) = ( cff1 * p_cosx2(i1) - cff2 ) * p_pmj1(i1) - p_pmj(i1)
      end do
    end do
    
    do i2 = n64+1, this%n, 16
      p_cosx2(1:16) => this%cosx2(i2:i2+15)
      p_pmj(1:16)   => pmj(i2:i2+15)
      p_pmj1(1:16)  => pmj1(i2:i2+15)
      
      !$omp simd aligned (p_pmj,p_pmj1,p_cosx2:alig)
      do i1 = 1, 16
        p_pmj(i1) = ( cff1 * p_cosx2(i1) - cff2 ) * p_pmj1(i1) - p_pmj(i1)
      end do
    end do
    
  end procedure mj_rec_sub
  
end submodule poly_mj