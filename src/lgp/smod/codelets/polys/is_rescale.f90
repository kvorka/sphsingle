submodule (lege_poly) is_rescale
  implicit none; contains
  
  module procedure is_rescale_sub
    integer        :: i1, i2
    real(kind=dbl) :: cff1, cff2, cff3, cff4
    
    do i2 = 1, (this%nrma/4)*4, 4
      cff1 = this%amj(i2  )
      cff2 = this%amj(i2+1)
      cff3 = this%amj(i2+2)
      cff4 = this%amj(i2+3)
      
      !$omp simd
      do i1 = 1, 4
        rcab(i1,i2  ) = cff1 * rcab(i1,i2  )
        rcab(i1,i2+1) = cff2 * rcab(i1,i2+1)
        rcab(i1,i2+2) = cff3 * rcab(i1,i2+2)
        rcab(i1,i2+3) = cff4 * rcab(i1,i2+3)
      end do
    end do
    
    do i2 = (this%nrma/4)*4+1, this%nrma
      cff1 = this%amj(i2)
      
      !$omp simd
      do i1 = 1, 4
        rcab(i1,i2) = cff1 * rcab(i1,i2)
      end do
    end do
    
  end procedure is_rescale_sub
  
end submodule is_rescale