submodule (lege_poly) bwd_sum
  implicit none; contains
  
  module procedure bwd_sum_sub
    integer        :: i1, i2
    real(kind=dbl) :: p1, p2, c1, c2, c3, c4
    
    c1 = cc(1)
    c2 = cc(2)
    c3 = cc(3)
    c4 = cc(4)
    
    do i2 = 1, this%n_dbl_2
      !$omp simd
      do i1 = 1, ndbl
        p1 = pmj(i1,1,i2)
        p2 = pmj(i1,2,i2)
        
        swork(i1,1,i2) = swork(i1,1,i2) + p1 * c1
        swork(i1,5,i2) = swork(i1,5,i2) + p2 * c1
        swork(i1,2,i2) = swork(i1,2,i2) + p1 * c2
        swork(i1,6,i2) = swork(i1,6,i2) + p2 * c2
        swork(i1,3,i2) = swork(i1,3,i2) + p1 * c3
        swork(i1,7,i2) = swork(i1,7,i2) + p2 * c3
        swork(i1,4,i2) = swork(i1,4,i2) + p1 * c4
        swork(i1,8,i2) = swork(i1,8,i2) + p2 * c4
      end do
    end do
    
  end procedure bwd_sum_sub
  
end submodule bwd_sum