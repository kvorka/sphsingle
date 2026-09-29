submodule (lege_poly) fwd_sum
  implicit none; contains
  
  module procedure fwd_sum_m_sub
    integer                             :: i1, i2, ima
    real(kind=dbl)                      :: f1, f2, cr1, cr2, cr3, cr4, p0, p1, p2, p3, s0, s1, s2, s3
    real(kind=dbl), pointer, contiguous :: pmj2(:,:,:)
    
    !!! After the FFT, we need to shuffle the packing north/south and real/imaginary
    !!! into packing suitable for summation.
    do i2 = 1, n1
      !$omp simd
      do i1 = 1, ndbl
        f1 = wght(i1,1,i2)
        f2 = cosx(i1,1,i2) * f1
        
        cr1 = grid(i1,1,i2,1)
        cr2 = grid(i1,1,i2,2)
        cr3 = grid(i1,1,i2,3)
        cr4 = grid(i1,1,i2,4)
        
        s0 = ( cr1 - cr2 ) * f1
        s2 = ( cr1 + cr2 ) * f2
        s1 = ( cr3 - cr4 ) * f1
        s3 = ( cr3 + cr4 ) * f2
        
        swork(i1,1,i2) = s0
        swork(i1,2,i2) = s1
        swork(i1,3,i2) = s2
        swork(i1,4,i2) = s3
        
        f1 = wght(i1,2,i2)
        f2 = cosx(i1,2,i2) * f1
        
        cr1 = grid(i1,2,i2,1)
        cr2 = grid(i1,2,i2,2)
        cr3 = grid(i1,2,i2,3)
        cr4 = grid(i1,2,i2,4)
        
        s0 = ( cr1 - cr2 ) * f1
        s2 = ( cr1 + cr2 ) * f2
        s1 = ( cr3 - cr4 ) * f1
        s3 = ( cr3 + cr4 ) * f2
        
        swork(i1,5,i2) = s0
        swork(i1,6,i2) = s1
        swork(i1,7,i2) = s2
        swork(i1,8,i2) = s3
        
        f1 = wght(i1,3,i2)
        f2 = cosx(i1,3,i2) * f1
        
        cr1 = grid(i1,3,i2,1)
        cr2 = grid(i1,3,i2,2)
        cr3 = grid(i1,3,i2,3)
        cr4 = grid(i1,3,i2,4)
        
        s0 = ( cr1 - cr2 ) * f1
        s2 = ( cr1 + cr2 ) * f2
        s1 = ( cr3 - cr4 ) * f1
        s3 = ( cr3 + cr4 ) * f2
        
        swork(i1, 9,i2) = s0
        swork(i1,10,i2) = s1
        swork(i1,11,i2) = s2
        swork(i1,12,i2) = s3
        
        f1 = wght(i1,4,i2)
        f2 = cosx(i1,4,i2) * f1
        
        cr1 = grid(i1,4,i2,1)
        cr2 = grid(i1,4,i2,2)
        cr3 = grid(i1,4,i2,3)
        cr4 = grid(i1,4,i2,4)
        
        s0 = ( cr1 - cr2 ) * f1
        s2 = ( cr1 + cr2 ) * f2
        s1 = ( cr3 - cr4 ) * f1
        s3 = ( cr3 + cr4 ) * f2
        
        swork(i1,13,i2) = s0
        swork(i1,14,i2) = s1
        swork(i1,15,i2) = s2
        swork(i1,16,i2) = s3
      end do
    end do
    
    !! Starting from degree j equal to order m, we need to forward the recursion for
    !! pmm, prepare the recursion for pmj by setting pmj1 to zero, and we need to
    !! set the initial value of swork to cc * pmj (first member of the sum).
    cr1 = zero
    cr2 = zero
    cr3 = zero
    cr4 = zero
    
    !$omp simd collapse (2) reduction (+:cr1,cr2,cr3,cr4)
    do i2 = 1, n1
      do i1 = 1, ndbl
        pmj1(i1,1,i2) = zero
        pmj1(i1,2,i2) = zero
        pmj1(i1,3,i2) = zero
        pmj1(i1,4,i2) = zero
        
        p0 = pmm(i1,1,i2)
        p1 = pmm(i1,2,i2)
        p2 = pmm(i1,3,i2)
        p3 = pmm(i1,4,i2)
        
        pmj(i1,1,i2) = p0
        pmj(i1,2,i2) = p1
        pmj(i1,3,i2) = p2
        pmj(i1,4,i2) = p3
        
        s0 = swork(i1,1,i2)
        s1 = swork(i1,2,i2)
        s2 = swork(i1,3,i2)
        s3 = swork(i1,4,i2)
        
        cr1 = cr1 + p0 * s0
        cr2 = cr2 + p0 * s1
        cr3 = cr3 + p0 * s2
        cr4 = cr4 + p0 * s3
        
        s0 = swork(i1,5,i2)
        s1 = swork(i1,6,i2)
        s2 = swork(i1,7,i2)
        s3 = swork(i1,8,i2)
        
        cr1 = cr1 + p1 * s0
        cr2 = cr2 + p1 * s1
        cr3 = cr3 + p1 * s2
        cr4 = cr4 + p1 * s3
        
        s0 = swork(i1, 9,i2)
        s1 = swork(i1,10,i2)
        s2 = swork(i1,11,i2)
        s3 = swork(i1,12,i2)
        
        cr1 = cr1 + p2 * s0
        cr2 = cr2 + p2 * s1
        cr3 = cr3 + p2 * s2
        cr4 = cr4 + p2 * s3
        
        s0 = swork(i1,13,i2)
        s1 = swork(i1,14,i2)
        s2 = swork(i1,15,i2)
        s3 = swork(i1,16,i2)
        
        cr1 = cr1 + p3 * s0
        cr2 = cr2 + p3 * s1
        cr3 = cr3 + p3 * s2
        cr4 = cr4 + p3 * s3
      end do
    end do
    
    cr(1,ma1) = cr(1,ma1) + cr1
    cr(2,ma1) = cr(2,ma1) + cr2
    cr(3,ma1) = cr(3,ma1) + cr3
    cr(4,ma1) = cr(4,ma1) + cr4
    
    !! Following with the recursion for degrees m+1 to jmax. We need to repointer our
    !! polynomials, follow with recursion and add cc * pmj to our swork accumulator.
    do ima = ma1+1, ma2
      f1 = fmj(1,ima)
      f2 = fmj(2,ima)
      
      cr1 = zero
      cr2 = zero
      cr3 = zero
      cr4 = zero
      
      pmj2 => pmj1
      pmj1 => pmj
      pmj  => pmj2
      
      !$omp simd collapse (2) reduction (+:cr1,cr2,cr3,cr4)
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
          
          s0 = cosx2(i1,3,i2)
          s1 = cosx2(i1,4,i2)
          s2 = pmj1(i1,3,i2)
          s3 = pmj1(i1,4,i2)
          
          s0 = f1 * s0 - f2
          s1 = f1 * s1 - f2
          
          p2 = pmj(i1,3,i2)
          p3 = pmj(i1,4,i2)
          
          p2 = s0 * s2 - p2
          p3 = s1 * s3 - p3
          
          pmj(i1,3,i2) = p2
          pmj(i1,4,i2) = p3
          
          s0 = swork(i1,1,i2)
          s1 = swork(i1,2,i2)
          s2 = swork(i1,3,i2)
          s3 = swork(i1,4,i2)
          
          cr1 = cr1 + p0 * s0
          cr2 = cr2 + p0 * s1
          cr3 = cr3 + p0 * s2
          cr4 = cr4 + p0 * s3
          
          s0 = swork(i1,5,i2)
          s1 = swork(i1,6,i2)
          s2 = swork(i1,7,i2)
          s3 = swork(i1,8,i2)
          
          cr1 = cr1 + p1 * s0
          cr2 = cr2 + p1 * s1
          cr3 = cr3 + p1 * s2
          cr4 = cr4 + p1 * s3
          
          s0 = swork(i1, 9,i2)
          s1 = swork(i1,10,i2)
          s2 = swork(i1,11,i2)
          s3 = swork(i1,12,i2)
          
          cr1 = cr1 + p2 * s0
          cr2 = cr2 + p2 * s1
          cr3 = cr3 + p2 * s2
          cr4 = cr4 + p2 * s3
          
          s0 = swork(i1,13,i2)
          s1 = swork(i1,14,i2)
          s2 = swork(i1,15,i2)
          s3 = swork(i1,16,i2)
          
          cr1 = cr1 + p3 * s0
          cr2 = cr2 + p3 * s1
          cr3 = cr3 + p3 * s2
          cr4 = cr4 + p3 * s3
        end do
      end do
      
      cr(1,ima) = cr(1,ima) + cr1
      cr(2,ima) = cr(2,ima) + cr2
      cr(3,ima) = cr(3,ima) + cr3
      cr(4,ima) = cr(4,ima) + cr4
    end do
    
  end procedure fwd_sum_m_sub
  
end submodule fwd_sum