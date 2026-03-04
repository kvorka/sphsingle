submodule (lege_poly) bwd_sum_m
  implicit none; contains
  
  module procedure bwd_sum_m_sub
    integer                             :: n2, i1, i2, ima
    real(kind=dbl)                      :: c1, c2, c3, c4, cff1, cff2
    real(kind=dbl), pointer, contiguous :: pmj2(:,:)
    
    !! As the cycles are going to be unrolled by hand by a factor of two for higher
    !! efficiency, we need n1/2 (will be used later to test, wheter n1 is divisible
    !! by two without a change).
    n2 = ( n1 / 2 ) * 2
    
    !! Starting from degree j equal to order m, we need to forward the recursion for
    !! pmm, prepare the recursion for pmj by setting pmj1 to zero, and we need to
    !! set the initial value of swork to cc * pmj (first member of the sum).
    c1 = cc(1,ma1)
    c2 = cc(2,ma1)
    c3 = cc(3,ma1)
    c4 = cc(4,ma1)
    
    cff2 = fmj(2,ma1)
    
    select case (ma1)
      case (1)
        do i2 = 1, n1
          !$omp simd aligned (pmj1:alig)
          do i1 = 1, ndbl
            pmm(i1,i2)  = cff2
            pmj1(i1,i2) = zero
          end do
        end do
        
      case default
        do i2 = 1, n1
          !$omp simd aligned (pmj1:alig)
          do i1 = 1, ndbl
            pmm(i1,i2)  = cff2 * sinx(i1,i2) * pmm(i1,i2)
            pmj1(i1,i2) = zero
          end do
        end do
    end select
    
    do i2 = 1, n2, 2
      !$omp simd aligned (pmj:alig)
      do i1 = 1, ndbl
        pmj(i1,i2  ) = pmm(i1,i2  ) / cosx(i1,i2  )
        pmj(i1,i2+1) = pmm(i1,i2+1) / cosx(i1,i2+1)
        
        swork(i1,1,i2) = pmj(i1,i2) * c1
        swork(i1,2,i2) = pmj(i1,i2) * c2
        swork(i1,3,i2) = pmj(i1,i2) * c3
        swork(i1,4,i2) = pmj(i1,i2) * c4
        
        swork(i1,1,i2+1) = pmj(i1,i2+1) * c1
        swork(i1,2,i2+1) = pmj(i1,i2+1) * c2
        swork(i1,3,i2+1) = pmj(i1,i2+1) * c3
        swork(i1,4,i2+1) = pmj(i1,i2+1) * c4
      end do
    end do
    
    if ( n2 /= n1 ) then
      !$omp simd aligned (pmj:alig)
      do i1 = 1, ndbl
        pmj(i1,n1) = pmm(i1,n1) / cosx(i1,n1)
        
        swork(i1,1,n1) = pmj(i1,n1) * c1
        swork(i1,2,n1) = pmj(i1,n1) * c2
        swork(i1,3,n1) = pmj(i1,n1) * c3
        swork(i1,4,n1) = pmj(i1,n1) * c4
      end do
    end if
    
    !! Following with the recursion for degrees m+1 to jmax. We need to repointer our
    !! polynomials, follow with recursion and add cc * pmj to our swork accumulator.
    do ima = ma1+1, ma2
      pmj2 => pmj1
      pmj1 => pmj
      pmj  => pmj2
      
      c1 = cc(1,ima)
      c2 = cc(2,ima)
      c3 = cc(3,ima)
      c4 = cc(4,ima)
      
      cff1 = fmj(1,ima)
      cff2 = fmj(2,ima)
      
      do i2 = 1, n2, 2
        !$omp simd aligned (pmj,pmj1:alig)
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
        !$omp simd aligned (pmj,pmj1:alig)
        do i1 = 1, ndbl
          pmj(i1,n1) = ( cff1 * cosx2(i1,n1) - cff2 ) * pmj1(i1,n1) - pmj(i1,n1)
          
          swork(i1,1,n1) = swork(i1,1,n1) + pmj(i1,n1) * c1
          swork(i1,2,n1) = swork(i1,2,n1) + pmj(i1,n1) * c2
          swork(i1,3,n1) = swork(i1,3,n1) + pmj(i1,n1) * c3
          swork(i1,4,n1) = swork(i1,4,n1) + pmj(i1,n1) * c4
        end do
      end if
    end do
    
    !! As we are done with computing the summation, we need to reshufle the data from
    !! packed sum to south/north and real/imaginary parts for upcomming FFT.
    do i2 = 1, n1
      !$omp simd
      do i1 = 1, ndbl
        cff1 = cosx(i1,i2)
        
        c1 = swork(i1,1,i2)
        c2 = swork(i1,3,i2)
        c3 = swork(i1,2,i2)
        c4 = swork(i1,4,i2)
        
        grid(i1,i2,1) = c2 * cff1 + c1
        grid(i1,i2,2) = c2 * cff1 - c1
        grid(i1,i2,3) = c4 * cff1 + c3
        grid(i1,i2,4) = c4 * cff1 - c3
      end do
    end do
    
  end procedure bwd_sum_m_sub
  
end submodule bwd_sum_m