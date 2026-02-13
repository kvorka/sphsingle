submodule (lege_poly) bwd_sum
  implicit none; contains
  
  module procedure bwd_sum_sub
    integer        :: i1, i2
    real(kind=dbl) :: c1, c2, c3, c4
    
    c1 = cc(1)
    c2 = cc(2)
    c3 = cc(3)
    c4 = cc(4)
    
    do i2 = 1, this%n_dbl, 2
      !$omp simd
      do i1 = 1, ndbl
        swork(i1,i2  ,1) = swork(i1,i2  ,1) + pmj(i1,i2  ) * c1
        swork(i1,i2+1,1) = swork(i1,i2+1,1) + pmj(i1,i2+1) * c1
        swork(i1,i2  ,2) = swork(i1,i2  ,2) + pmj(i1,i2  ) * c2
        swork(i1,i2+1,2) = swork(i1,i2+1,2) + pmj(i1,i2+1) * c2
        swork(i1,i2  ,3) = swork(i1,i2  ,3) + pmj(i1,i2  ) * c3
        swork(i1,i2+1,3) = swork(i1,i2+1,3) + pmj(i1,i2+1) * c3
        swork(i1,i2  ,4) = swork(i1,i2  ,4) + pmj(i1,i2  ) * c4
        swork(i1,i2+1,4) = swork(i1,i2+1,4) + pmj(i1,i2+1) * c4
      end do
    end do
    
  end procedure bwd_sum_sub
  
end submodule bwd_sum