submodule (lege_poly) fwd_sum2
  implicit none; contains
  
  module procedure fwd_sum2_sub
    integer        :: i1, i2, n2
    real(kind=dbl) :: cff1, cff2
    
    cff1 = fmj(1)
    cff2 = fmj(2)
    
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
    
    n2 = ( n1 / 2 ) * 2
    
    do i2 = 1, n2, 2
      !$omp simd
      do i1 = 1, ndbl
        pmj(i1,i2  ) = ( cff1 * cosx2(i1,i2  ) - cff2 ) * pmj1(i1,i2  ) - pmj(i1,i2  )
        pmj(i1,i2+1) = ( cff1 * cosx2(i1,i2+1) - cff2 ) * pmj1(i1,i2+1) - pmj(i1,i2+1)
        
        acc(i1,1) = acc(i1,1) + pmj(i1,i2) * s1(i1,i2)
        acc(i1,2) = acc(i1,2) + pmj(i1,i2) * s2(i1,i2)
        acc(i1,3) = acc(i1,3) + pmj(i1,i2) * s3(i1,i2)
        acc(i1,4) = acc(i1,4) + pmj(i1,i2) * s4(i1,i2)
        
        acc2(i1,1) = acc2(i1,1) + pmj(i1,i2+1) * s1(i1,i2+1)
        acc2(i1,2) = acc2(i1,2) + pmj(i1,i2+1) * s2(i1,i2+1)
        acc2(i1,3) = acc2(i1,3) + pmj(i1,i2+1) * s3(i1,i2+1)
        acc2(i1,4) = acc2(i1,4) + pmj(i1,i2+1) * s4(i1,i2+1)
      end do
    end do
    
    if ( n2 /= n1 ) then
      !$omp simd
      do i1 = 1, ndbl
        pmj(i1,n1) = ( cff1 * cosx2(i1,n1) - cff2 ) * pmj1(i1,n1) - pmj(i1,n1)
        
        acc(i1,1) = acc(i1,1) + pmj(i1,n1) * s1(i1,n1)
        acc(i1,2) = acc(i1,2) + pmj(i1,n1) * s2(i1,n1)
        acc(i1,3) = acc(i1,3) + pmj(i1,n1) * s3(i1,n1)
        acc(i1,4) = acc(i1,4) + pmj(i1,n1) * s4(i1,n1)
      end do
    end if
    
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
    
  end procedure fwd_sum2_sub
  
end submodule fwd_sum2