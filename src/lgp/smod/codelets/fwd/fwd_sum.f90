submodule (lege_poly) fwd_sum
  implicit none; contains
  
  module procedure fwd_sum_sub
    integer        :: i1, i2
    real(kind=dbl) :: p1, p2
    
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
      !$omp simd
      do i1 = 1, ndbl
        p1 = pmj(i1,1,i2)
        p2 = pmj(i1,2,i2)
        
        acc(i1,1) = acc(i1,1) + p1 * swork(i1,1,i2,1)
        acc(i1,2) = acc(i1,2) + p1 * swork(i1,1,i2,2)
        acc(i1,3) = acc(i1,3) + p1 * swork(i1,1,i2,3)
        acc(i1,4) = acc(i1,4) + p1 * swork(i1,1,i2,4)
        
        acc2(i1,1) = acc2(i1,1) + p2 * swork(i1,2,i2,1)
        acc2(i1,2) = acc2(i1,2) + p2 * swork(i1,2,i2,2)
        acc2(i1,3) = acc2(i1,3) + p2 * swork(i1,2,i2,3)
        acc2(i1,4) = acc2(i1,4) + p2 * swork(i1,2,i2,4)
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
    
  end procedure fwd_sum_sub
  
end submodule fwd_sum
