submodule (lege_poly) bwd_sum
  implicit none; contains
  
  module procedure bwd_sum_m_sub
    integer                             :: i1, i2, ima
    real(kind=dbl)                      :: f1, f2, c1, c2, c3, c4, p0, p1, s0, s1, s2, s3
    real(kind=dbl), pointer, contiguous :: pmj2(:,:,:)
    
    !! Starting from degree j equal to order m, we need to forward the recursion for
    !! pmm, prepare the recursion for pmj by setting pmj1 to zero, and we need to
    !! set the initial value of swork to cc * pmj (first member of the sum).
    c1 = cc(1,ma1)
    c2 = cc(2,ma1)
    c3 = cc(3,ma1)
    c4 = cc(4,ma1)
    
    !$omp simd collapse (2)
    do i2 = 1, n1
      do i1 = 1, ndbl
        pmj1(i1,1,i2) = zero
        pmj1(i1,2,i2) = zero
        pmj1(i1,3,i2) = zero
        pmj1(i1,4,i2) = zero
        
        p0 = pmm(i1,1,i2)
        p1 = pmm(i1,2,i2)
        
        pmj(i1,1,i2) = p0
        pmj(i1,2,i2) = p1
        
        s0 = p0 * c1
        s1 = p0 * c2
        s2 = p0 * c3
        s3 = p0 * c4
        
        swork(i1,1,i2) = s0
        swork(i1,2,i2) = s1
        swork(i1,3,i2) = s2
        swork(i1,4,i2) = s3
        
        s0 = p1 * c1
        s1 = p1 * c2
        s2 = p1 * c3
        s3 = p1 * c4
        
        swork(i1,5,i2) = s0
        swork(i1,6,i2) = s1
        swork(i1,7,i2) = s2
        swork(i1,8,i2) = s3
        
        p0 = pmm(i1,3,i2)
        p1 = pmm(i1,4,i2)
        
        pmj(i1,3,i2) = p0
        pmj(i1,4,i2) = p1
        
        s0 = p0 * c1
        s1 = p0 * c2
        s2 = p0 * c3
        s3 = p0 * c4
        
        swork(i1, 9,i2) = s0
        swork(i1,10,i2) = s1
        swork(i1,11,i2) = s2
        swork(i1,12,i2) = s3
        
        s0 = p1 * c1
        s1 = p1 * c2
        s2 = p1 * c3
        s3 = p1 * c4
        
        swork(i1,13,i2) = s0
        swork(i1,14,i2) = s1
        swork(i1,15,i2) = s2
        swork(i1,16,i2) = s3
      end do
    end do
    
    !! Following with the recursion for degrees m+1 to jmax. We need to repointer our
    !! polynomials, follow with recursion and add cc * pmj to our swork accumulator.
    do ima = ma1+1, ma2
      c1 = cc(1,ima)
      c2 = cc(2,ima)
      c3 = cc(3,ima)
      c4 = cc(4,ima)
      
      f1 = fmj(1,ima)
      f2 = fmj(2,ima)
      
      pmj2 => pmj1
      pmj1 => pmj
      pmj  => pmj2
      
      !$omp simd collapse (2)
      do i2 = 1, n1
        do i1 = 1, ndbl
          s0 = cosx2(i1,1,i2)
          s1 = cosx2(i1,2,i2)
          s2 = pmj1(i1,1,i2)
          s3 = pmj1(i1,2,i2)
          
          s0 = f1 * s0 - f2
          s1 = f1 * s1 - f2
          
          p0 = pmj(i1,1,i2)
          p1 = pmj(i1,2,i2)
          
          p0 = s0 * s2 - p0
          p1 = s1 * s3 - p1
          
          pmj(i1,1,i2) = p0
          pmj(i1,2,i2) = p1
          
          s0 = swork(i1,1,i2)
          s1 = swork(i1,2,i2)
          s2 = swork(i1,3,i2)
          s3 = swork(i1,4,i2)
          
          s0 = s0 + p0 * c1
          s1 = s1 + p0 * c2
          s2 = s2 + p0 * c3
          s3 = s3 + p0 * c4
          
          swork(i1,1,i2) = s0
          swork(i1,2,i2) = s1
          swork(i1,3,i2) = s2
          swork(i1,4,i2) = s3
          
          s0 = swork(i1,5,i2)
          s1 = swork(i1,6,i2)
          s2 = swork(i1,7,i2)
          s3 = swork(i1,8,i2)
          
          s0 = s0 + p1 * c1
          s1 = s1 + p1 * c2
          s2 = s2 + p1 * c3
          s3 = s3 + p1 * c4
          
          swork(i1,5,i2) = s0
          swork(i1,6,i2) = s1
          swork(i1,7,i2) = s2
          swork(i1,8,i2) = s3
          
          s0 = cosx2(i1,3,i2)
          s1 = cosx2(i1,4,i2)
          s2 = pmj1(i1,3,i2)
          s3 = pmj1(i1,4,i2)
          
          s0 = f1 * s0 - f2
          s1 = f1 * s1 - f2
          
          p0 = pmj(i1,3,i2)
          p1 = pmj(i1,4,i2)
          
          p0 = s0 * s2 - p0
          p1 = s1 * s3 - p1
          
          pmj(i1,3,i2) = p0
          pmj(i1,4,i2) = p1
          
          s0 = swork(i1, 9,i2)
          s1 = swork(i1,10,i2)
          s2 = swork(i1,11,i2)
          s3 = swork(i1,12,i2)
          
          s0 = s0 + p0 * c1
          s1 = s1 + p0 * c2
          s2 = s2 + p0 * c3
          s3 = s3 + p0 * c4
          
          swork(i1, 9,i2) = s0
          swork(i1,10,i2) = s1
          swork(i1,11,i2) = s2
          swork(i1,12,i2) = s3
          
          s0 = swork(i1,13,i2)
          s1 = swork(i1,14,i2)
          s2 = swork(i1,15,i2)
          s3 = swork(i1,16,i2)
          
          s0 = s0 + p1 * c1
          s1 = s1 + p1 * c2
          s2 = s2 + p1 * c3
          s3 = s3 + p1 * c4
          
          swork(i1,13,i2) = s0
          swork(i1,14,i2) = s1
          swork(i1,15,i2) = s2
          swork(i1,16,i2) = s3
        end do
      end do
    end do
    
    !! As we are done with computing the summation, we need to reshufle the data from
    !! packed sum to south/north and real/imaginary parts for upcomming FFT..
    !$omp simd collapse (2)
    do i2 = 1, n1
      do i1 = 1, ndbl
        p0 = cosx(i1,1,i2)
        p1 = cosx(i1,2,i2)
        
        s0 = swork(i1,1,i2)
        s1 = swork(i1,2,i2)
        s2 = swork(i1,3,i2)
        s3 = swork(i1,4,i2)
        
        c1 = s2 * p0 + s0
        c2 = s2 * p0 - s0
        c3 = s3 * p0 + s1
        c4 = s3 * p0 - s1
        
        grid(i1,1,i2,1) = c1
        grid(i1,1,i2,2) = c2
        grid(i1,1,i2,3) = c3
        grid(i1,1,i2,4) = c4
        
        s0 = swork(i1,5,i2)
        s1 = swork(i1,6,i2)
        s2 = swork(i1,7,i2)
        s3 = swork(i1,8,i2)
        
        c1 = s2 * p1 + s0
        c2 = s2 * p1 - s0
        c3 = s3 * p1 + s1
        c4 = s3 * p1 - s1
        
        grid(i1,2,i2,1) = c1
        grid(i1,2,i2,2) = c2
        grid(i1,2,i2,3) = c3
        grid(i1,2,i2,4) = c4
        
        p0 = cosx(i1,3,i2)
        p1 = cosx(i1,4,i2)
        
        s0 = swork(i1, 9,i2)
        s1 = swork(i1,10,i2)
        s2 = swork(i1,11,i2)
        s3 = swork(i1,12,i2)
        
        c1 = s2 * p0 + s0
        c2 = s2 * p0 - s0
        c3 = s3 * p0 + s1
        c4 = s3 * p0 - s1
        
        grid(i1,3,i2,1) = c1
        grid(i1,3,i2,2) = c2
        grid(i1,3,i2,3) = c3
        grid(i1,3,i2,4) = c4
        
        s0 = swork(i1,13,i2)
        s1 = swork(i1,14,i2)
        s2 = swork(i1,15,i2)
        s3 = swork(i1,16,i2)
        
        c1 = s2 * p1 + s0
        c2 = s2 * p1 - s0
        c3 = s3 * p1 + s1
        c4 = s3 * p1 - s1
        
        grid(i1,4,i2,1) = c1
        grid(i1,4,i2,2) = c2
        grid(i1,4,i2,3) = c3
        grid(i1,4,i2,4) = c4
      end do
    end do
    
  end procedure bwd_sum_m_sub
  
end submodule bwd_sum