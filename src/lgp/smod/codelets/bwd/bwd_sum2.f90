submodule (lege_poly) bwd_sum2
  implicit none; contains
  
  module procedure bwd_sum2_sub
    integer        :: i1, i2, n2
    real(kind=dbl) :: c1, c2, c3, c4, cff1, cff2
    
    n2 = ( n1 / 2 ) * 2
    
    c1 = cc(1)
    c2 = cc(2)
    c3 = cc(3)
    c4 = cc(4)
    
    cff1 = fmj(1)
    cff2 = fmj(2)
    
    do i2 = 1, n2, 2
      !$omp simd
      do i1 = 1, ndbl
        pmj(i1,i2  ) = ( cff1 * cosx2(i1,i2  ) - cff2 ) * pmj1(i1,i2  ) - pmj(i1,i2  )
        pmj(i1,i2+1) = ( cff1 * cosx2(i1,i2+1) - cff2 ) * pmj1(i1,i2+1) - pmj(i1,i2+1)
        
        swork(i1,1,i2) = swork(i1,1,i2) + pmj(i1,i2) * c1
        swork(i1,2,i2) = swork(i1,2,i2) + pmj(i1,i2) * c2
        swork(i1,3,i2) = swork(i1,3,i2) + pmj(i1,i2) * c3
        swork(i1,4,i2) = swork(i1,4,i2) + pmj(i1,i2) * c4
        
        swork(i1,1,i2+1) = swork(i1,1,i2+1) + pmj(i1,i2+1) * c1
        swork(i1,2,i2+1) = swork(i1,2,i2+1) + pmj(i1,i2+1) * c2
        swork(i1,3,i2+1) = swork(i1,3,i2+1) + pmj(i1,i2+1) * c3
        swork(i1,4,i2+1) = swork(i1,4,i2+1) + pmj(i1,i2+1) * c4
      end do
    end do
    
    if ( n2 /= n1 ) then
      !$omp simd
      do i1 = 1, ndbl
        pmj(i1,n1) = ( cff1 * cosx2(i1,n1) - cff2 ) * pmj1(i1,n1) - pmj(i1,n1)
        
        swork(i1,1,n1) = swork(i1,1,n1) + pmj(i1,n1) * c1
        swork(i1,2,n1) = swork(i1,2,n1) + pmj(i1,n1) * c2
        swork(i1,3,n1) = swork(i1,3,n1) + pmj(i1,n1) * c3
        swork(i1,4,n1) = swork(i1,4,n1) + pmj(i1,n1) * c4
      end do
    end if
    
  end procedure bwd_sum2_sub
  
end submodule bwd_sum2