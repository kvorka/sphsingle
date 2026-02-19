submodule (lege_poly) bwd_sum1
  implicit none; contains
  
  module procedure bwd_sum1_sub
    integer        :: i1, i2, n2
    real(kind=dbl) :: c1, c2, c3, c4
    
    select case (ma)
      case (1)
        do i2 = 1, n1
          !$omp simd
          do i1 = 1, ndbl
            pmm(i1,i2)  = fmj
            pmj1(i1,i2) = zero
          end do
        end do
        
      case default
        do i2 = 1, n1
          !$omp simd
          do i1 = 1, ndbl
            pmm(i1,i2)  = fmj * sinx(i1,i2) * pmm(i1,i2)
            pmj1(i1,i2) = zero
          end do
        end do
    end select
    
    c1 = cc(1)
    c2 = cc(2)
    c3 = cc(3)
    c4 = cc(4)
    
    n2 = ( n1 / 2 ) * 2
    
    do i2 = 1, n2, 2
      !$omp simd
      do i1 = 1, ndbl
        pmj(i1,i2  ) = pmm(i1,i2  ) / cosx(i1,i2  )
        pmj(i1,i2+1) = pmm(i1,i2+1) / cosx(i1,i2+1)
        
        s1(i1,i2) = pmj(i1,i2) * c1
        s2(i1,i2) = pmj(i1,i2) * c2
        s3(i1,i2) = pmj(i1,i2) * c3
        s4(i1,i2) = pmj(i1,i2) * c4
        
        s1(i1,i2+1) = pmj(i1,i2+1) * c1
        s2(i1,i2+1) = pmj(i1,i2+1) * c2
        s3(i1,i2+1) = pmj(i1,i2+1) * c3
        s4(i1,i2+1) = pmj(i1,i2+1) * c4
      end do
    end do
    
    if ( n2 /= n1 ) then
      !$omp simd
      do i1 = 1, ndbl
        pmj(i1,n1) = pmm(i1,n1) / cosx(i1,n1)
        
        s1(i1,n1) = pmj(i1,n1) * c1
        s2(i1,n1) = pmj(i1,n1) * c2
        s3(i1,n1) = pmj(i1,n1) * c3
        s4(i1,n1) = pmj(i1,n1) * c4
      end do
    end if
    
  end procedure bwd_sum1_sub
  
end submodule bwd_sum1