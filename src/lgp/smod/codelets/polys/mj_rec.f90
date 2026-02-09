submodule (lege_poly) poly_mj
  implicit none; contains
  
  module procedure mj_rec_sub
    integer                             :: i1, i2
    real(kind=dbl)                      :: cff1, cff2, fac1, fac2, pmjA, pmjB
    real(kind=dbl), pointer, contiguous :: p1_pmj(:), p1_pmj1(:), p1_cosx2(:), &
                                         & p2_pmj(:), p2_pmj1(:), p2_cosx2(:)
    
    cff1 = this%fmj(1,ima)
    cff2 = this%fmj(2,ima)
    
    do i2 = 1, this%n32, 32
      p1_cosx2(1:16) => this%cosx2(i2   :i2+15)
      p2_cosx2(1:16) => this%cosx2(i2+16:i2+31)
      
      p1_pmj(1:16) => pmj(i2   :i2+15)
      p2_pmj(1:16) => pmj(i2+16:i2+31)
      
      p1_pmj1(1:16) => pmj1(i2   :i2+15)
      p2_pmj1(1:16) => pmj1(i2+16:i2+31)
      
      !$omp simd aligned (p1_pmj,p1_pmj1,p1_cosx2,p2_pmj,p2_pmj1,p2_cosx2:alig)
      do i1 = 1, 16
        fac1 = cff1 * p1_cosx2(i1) - cff2
        fac2 = cff1 * p2_cosx2(i1) - cff2
        
        pmjA = p1_pmj(i1)
        pmjB = p2_pmj(i1)
        
        p1_pmj(i1) = ( cff1 * p1_cosx2(i1) - cff2 ) * p1_pmj1(i1) - pmjA
        p2_pmj(i1) = ( cff1 * p2_cosx2(i1) - cff2 ) * p2_pmj1(i1) - pmjB
      end do
    end do
    
    do i2 = this%n32+1, this%n, 16
      p1_cosx2(1:16) => this%cosx2(i2:i2+15)
      p1_pmj(1:16)   => pmj(i2:i2+15)
      p1_pmj1(1:16)  => pmj1(i2:i2+15)
      
      !$omp simd aligned (p1_pmj,p1_pmj1,p1_cosx2:alig)
      do i1 = 1, 16
        p1_pmj(i1) = ( cff1 * p1_cosx2(i1) - cff2 ) * p1_pmj1(i1) - p1_pmj(i1)
      end do
    end do
    
  end procedure mj_rec_sub
  
end submodule poly_mj