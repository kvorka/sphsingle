submodule (lege_poly) fwd_sum
  implicit none; contains
  
  module procedure fwd_sum_m_sub
    integer                             :: i1, i2, ima, n2
    real(kind=dbl)                      :: cff1, cff2, g1, g2, g3, g4
    real(kind=dbl), pointer, contiguous :: pmj2(:,:)
    
    !! As the cycles are going to be unrolled by hand by a factor of two for higher
    !! efficiency, we need n1/2 (will be used later to test, wheter n1 is divisible
    !! by two without a change).
    n2 = ( n1 / 2 ) * 2
    
    !!! After the FFT, we need to shuffle the packing north/south and real/imaginary
    !!! into packing suitable for summation.
    !GCC$ unroll 4
    !DIR$ unroll (4)
    do i2 = 1, n1
      !$omp simd
      do i1 = 1, ndbl
        cff1 = wght(i1,i2)
        cff2 = cosx(i1,i2)
        
        g1 = grid(i1,i2,1)
        g2 = grid(i1,i2,2)
        g3 = grid(i1,i2,3)
        g4 = grid(i1,i2,4)
        
        swork(i1,1,i2) = ( g1 - g2 ) * cff1
        swork(i1,3,i2) = ( g1 + g2 ) * cff1 * cff2
        swork(i1,2,i2) = ( g3 - g4 ) * cff1
        swork(i1,4,i2) = ( g3 + g4 ) * cff1 * cff2
      end do
    end do
    
    !! Starting from degree j equal to order m, we need to forward the recursion for
    !! pmm, prepare the recursion for pmj by setting pmj1 to zero, and we need to
    !! set the initial value of swork to cc * pmj (first member of the sum).
    !$omp simd
    do i1 = 1, ndbl
      !GCC$ unroll 4
      !DIR$ unroll (4)
      do i2 = 1, 4
        acc(i1,i2)  = zero
        acc2(i1,i2) = zero
      end do
    end do
    
    do i2 = 1, n2, 2
      !$omp simd aligned (pmj1,pmj,pmm:alig)
      do i1 = 1, ndbl
        pmj1(i1,i2  ) = zero
        pmj1(i1,i2+1) = zero
        
        pmj(i1,i2  ) = pmm(i1,i2  )
        pmj(i1,i2+1) = pmm(i1,i2+1)
        
        acc(i1,1) = acc(i1,1) + pmj(i1,i2) * swork(i1,1,i2)
        acc(i1,2) = acc(i1,2) + pmj(i1,i2) * swork(i1,2,i2)
        acc(i1,3) = acc(i1,3) + pmj(i1,i2) * swork(i1,3,i2)
        acc(i1,4) = acc(i1,4) + pmj(i1,i2) * swork(i1,4,i2)
        
        acc2(i1,1) = acc2(i1,1) + pmj(i1,i2+1) * swork(i1,1,i2+1)
        acc2(i1,2) = acc2(i1,2) + pmj(i1,i2+1) * swork(i1,2,i2+1)
        acc2(i1,3) = acc2(i1,3) + pmj(i1,i2+1) * swork(i1,3,i2+1)
        acc2(i1,4) = acc2(i1,4) + pmj(i1,i2+1) * swork(i1,4,i2+1)
      end do
    end do
    
    if ( n2 /= n1 ) then
      !$omp simd aligned (pmj1,pmj,pmm:alig)
      do i1 = 1, ndbl
        pmj1(i1,n1) = zero
        
        pmj(i1,n1) = pmm(i1,n1)
        
        acc(i1,1) = acc(i1,1) + pmj(i1,n1) * swork(i1,1,n1)
        acc(i1,2) = acc(i1,2) + pmj(i1,n1) * swork(i1,2,n1)
        acc(i1,3) = acc(i1,3) + pmj(i1,n1) * swork(i1,3,n1)
        acc(i1,4) = acc(i1,4) + pmj(i1,n1) * swork(i1,4,n1)
      end do
    end if
    
    !$omp simd
    do i1 = 1, ndbl
      cr(1,ma1) = cr(1,ma1) + acc(i1,1) + acc2(i1,1)
      cr(2,ma1) = cr(2,ma1) + acc(i1,2) + acc2(i1,2)
      cr(3,ma1) = cr(3,ma1) + acc(i1,3) + acc2(i1,3)
      cr(4,ma1) = cr(4,ma1) + acc(i1,4) + acc2(i1,4)
    end do
    
    !! Following with the recursion for degrees m+1 to jmax. We need to repointer our
    !! polynomials, follow with recursion and add cc * pmj to our swork accumulator.
    do ima = ma1+1, ma2
      cff1 = fmj(1,ima)
      cff2 = fmj(2,ima)
      
      !$omp simd
      do i1 = 1, ndbl
        !GCC$ unroll 4
        !DIR$ unroll (4)
        do i2 = 1, 4
          acc(i1,i2)  = zero
          acc2(i1,i2) = zero
        end do
      end do
      
      pmj2 => pmj1
      pmj1 => pmj
      pmj  => pmj2
      
      do i2 = 1, n2, 2
        !$omp simd aligned (pmj1,pmj:alig)
        do i1 = 1, ndbl
          pmj(i1,i2  ) = ( cff1 * cosx2(i1,i2  ) - cff2 ) * pmj1(i1,i2  ) - pmj(i1,i2  )
          pmj(i1,i2+1) = ( cff1 * cosx2(i1,i2+1) - cff2 ) * pmj1(i1,i2+1) - pmj(i1,i2+1)
          
          acc(i1,1) = acc(i1,1) + pmj(i1,i2) * swork(i1,1,i2)
          acc(i1,2) = acc(i1,2) + pmj(i1,i2) * swork(i1,2,i2)
          acc(i1,3) = acc(i1,3) + pmj(i1,i2) * swork(i1,3,i2)
          acc(i1,4) = acc(i1,4) + pmj(i1,i2) * swork(i1,4,i2)
          
          acc2(i1,1) = acc2(i1,1) + pmj(i1,i2+1) * swork(i1,1,i2+1)
          acc2(i1,2) = acc2(i1,2) + pmj(i1,i2+1) * swork(i1,2,i2+1)
          acc2(i1,3) = acc2(i1,3) + pmj(i1,i2+1) * swork(i1,3,i2+1)
          acc2(i1,4) = acc2(i1,4) + pmj(i1,i2+1) * swork(i1,4,i2+1)
        end do
      end do
      
      if ( n2 /= n1 ) then
        !$omp simd aligned (pmj1,pmj:alig)
        do i1 = 1, ndbl
          pmj(i1,n1) = ( cff1 * cosx2(i1,n1) - cff2 ) * pmj1(i1,n1) - pmj(i1,n1)
          
          acc(i1,1) = acc(i1,1) + pmj(i1,n1) * swork(i1,1,n1)
          acc(i1,2) = acc(i1,2) + pmj(i1,n1) * swork(i1,2,n1)
          acc(i1,3) = acc(i1,3) + pmj(i1,n1) * swork(i1,3,n1)
          acc(i1,4) = acc(i1,4) + pmj(i1,n1) * swork(i1,4,n1)
        end do
      end if
      
      !$omp simd
      do i1 = 1, ndbl
        cr(1,ima) = cr(1,ima) + acc(i1,1) + acc2(i1,1)
        cr(2,ima) = cr(2,ima) + acc(i1,2) + acc2(i1,2)
        cr(3,ima) = cr(3,ima) + acc(i1,3) + acc2(i1,3)
        cr(4,ima) = cr(4,ima) + acc(i1,4) + acc2(i1,4)
      end do
    end do
    
  end procedure fwd_sum_m_sub
  
end submodule fwd_sum