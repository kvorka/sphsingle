submodule (lege_poly) bwd_sum
  implicit none; contains
  
  module procedure bwd_sum_sub
    integer        :: i1, i2
    real(kind=dbl) :: p1, p2, p3, p4, c1, c2, c3, c4
    
    c1 = cc(1)
    c2 = cc(2)
    c3 = cc(3)
    c4 = cc(4)
    
    
    do i2 = 1, this%n_dbl4, 4
      !$omp simd
      do i1 = 1, ndbl
        p1 = pmj(i1,i2  )
        p2 = pmj(i1,i2+1)
        p3 = pmj(i1,i2+2)
        p4 = pmj(i1,i2+3)

        swork(i1,i2  ,1) = swork(i1,i2  ,1) + p1 * c1
        swork(i1,i2+1,1) = swork(i1,i2+1,1) + p2 * c1
        swork(i1,i2+2,1) = swork(i1,i2+2,1) + p3 * c1
        swork(i1,i2+3,1) = swork(i1,i2+3,1) + p4 * c1
        
        swork(i1,i2  ,2) = swork(i1,i2  ,2) + p1 * c2
        swork(i1,i2+1,2) = swork(i1,i2+1,2) + p2 * c2
        swork(i1,i2+2,2) = swork(i1,i2+2,2) + p3 * c2
        swork(i1,i2+3,2) = swork(i1,i2+3,2) + p4 * c2

        swork(i1,i2  ,3) = swork(i1,i2  ,3) + p1 * c3
        swork(i1,i2+1,3) = swork(i1,i2+1,3) + p2 * c3
        swork(i1,i2+2,3) = swork(i1,i2+2,3) + p3 * c3
        swork(i1,i2+3,3) = swork(i1,i2+3,3) + p4 * c3
        
        swork(i1,i2  ,4) = swork(i1,i2  ,4) + p1 * c4
        swork(i1,i2+1,4) = swork(i1,i2+1,4) + p2 * c4
        swork(i1,i2+2,4) = swork(i1,i2+2,4) + p3 * c4
        swork(i1,i2+3,4) = swork(i1,i2+3,4) + p4 * c4
      end do
    end do
    
    if ( mod(this%n_dbl,4) /= 0 ) then
      i2 = this%n_dbl4+1
      
      !$omp simd
      do i1 = 1, ndbl
        p1 = pmj(i1,i2  )
        p2 = pmj(i1,i2+1)
        
        swork(i1,i2  ,1) = swork(i1,i2  ,1) + p1 * c1
        swork(i1,i2+1,1) = swork(i1,i2+1,1) + p2 * c1
        
        swork(i1,i2  ,2) = swork(i1,i2  ,2) + p1 * c2
        swork(i1,i2+1,2) = swork(i1,i2+1,2) + p2 * c2
        
        swork(i1,i2  ,3) = swork(i1,i2  ,3) + p1 * c3
        swork(i1,i2+1,3) = swork(i1,i2+1,3) + p2 * c3
        
        swork(i1,i2  ,4) = swork(i1,i2  ,4) + p1 * c4
        swork(i1,i2+1,4) = swork(i1,i2+1,4) + p2 * c4
      end do
    end if
    
  end procedure bwd_sum_sub
  
end submodule bwd_sum