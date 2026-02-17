submodule (lege_poly) fwd_sum1
  implicit none; contains
  
  module procedure fwd_sum1_sub
    integer :: i1, i2
    
    select case (ma)
      case (1)
        do i2 = 1, n1
          !$omp simd
          do i1 = 1, ndbl
            pmm(i1,1,i2) = fmj
            pmm(i1,2,i2) = fmj
            
            pmj1(i1,1,i2) = zero
            pmj1(i1,2,i2) = zero
          end do
        end do
        
      case default
        do i2 = 1, n1
          !$omp simd
          do i1 = 1, ndbl
            pmm(i1,1,i2)  = fmj * sinx(i1,1,i2) * pmm(i1,1,i2)
            pmm(i1,2,i2)  = fmj * sinx(i1,2,i2) * pmm(i1,2,i2)
            
            pmj1(i1,1,i2) = zero
            pmj1(i1,2,i2) = zero
          end do
        end do
    end select
    
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
    
    do i2 = 1, n1
      !$omp simd
      do i1 = 1, ndbl
        pmj(i1,1,i2) = pmm(i1,1,i2) / cosx(i1,1,i2)
        pmj(i1,2,i2) = pmm(i1,2,i2) / cosx(i1,2,i2)
        
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
    
  end procedure fwd_sum1_sub
  
end submodule fwd_sum1
