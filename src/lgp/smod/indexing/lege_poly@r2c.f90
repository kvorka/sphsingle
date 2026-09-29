submodule (lege_poly) r2c
  implicit none; contains
  
  module procedure index_fwd_sub
    integer                     :: m, j, jm, mj, ma
    real(kind=dbl), allocatable :: cab(:,:)
    
    !! First rescale required by the Ishioka recursion. This call is not
    !! inlined because the same rescale is required before the backword
    !! transform. Vectorized inside.
    call is_rescale_sub( this%nrma, this%amj, rcab )
    
    !! Reindexing after transform. Quadruplets of rescaled coefficients
    !! are synthethysed. These include real/imaginary, odd degree/even degree
    !! components, respectively.
    allocate( cab(2,this%jms) )
    
    m = 0
      !j == m
        ma = 1
        mj = 1
        
        cab(1,mj) = rcab(3,ma)
        cab(2,mj) = rcab(4,ma)
      
      !$omp simd
      do j = 1, (this%jmax-1)/2
        ma = ma+1
        mj = mj+2
        
        cab(1,mj-1) = this%emj(mj-1) * rcab(1,ma-1) + this%emj(mj) * rcab(1,ma)
        cab(2,mj-1) = this%emj(mj-1) * rcab(2,ma-1) + this%emj(mj) * rcab(2,ma)
        cab(1,mj  ) =                  rcab(3,ma  )
        cab(2,mj  ) =                  rcab(4,ma  )
      end do
      
      !j == jmax
      if ( mod(this%jmax,2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        cab(1,mj-1) = this%emj(mj-1) * rcab(1,ma-1) + this%emj(mj) * rcab(1,ma)
        cab(2,mj-1) = this%emj(mj-1) * rcab(2,ma-1) + this%emj(mj) * rcab(2,ma)
        cab(1,mj  ) =                  rcab(3,ma  )
        cab(2,mj  ) =                  rcab(4,ma  )
      
      else
        ma = ma+1
        mj = mj+1
        
        cab(1,mj) = this%emj(mj) * rcab(1,ma-1) + this%emj(mj+1) * rcab(1,ma)
        cab(2,mj) = this%emj(mj) * rcab(2,ma-1) + this%emj(mj+1) * rcab(2,ma)
      end if
    
    do m = 1, this%jmax-1
      !j == m
        ma = ma+1
        mj = mj+1
        
        cab(1,mj) = rcab(3,ma)
        cab(2,mj) = rcab(4,ma)
      
      !$omp simd
      do j = 1, (this%jmax-m-1)/2
        ma = ma+1
        mj = mj+2
        
        cab(1,mj-1) = this%emj(mj+m-1) * rcab(1,ma-1) + this%emj(mj+m) * rcab(1,ma)
        cab(2,mj-1) = this%emj(mj+m-1) * rcab(2,ma-1) + this%emj(mj+m) * rcab(2,ma)
        cab(1,mj  ) =                    rcab(3,ma  )
        cab(2,mj  ) =                    rcab(4,ma  )
      end do
      
      !j == jmax
      if ( mod((this%jmax-m),2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        cab(1,mj-1) = this%emj(mj+m-1) * rcab(1,ma-1) + this%emj(mj+m) * rcab(1,ma)
        cab(2,mj-1) = this%emj(mj+m-1) * rcab(2,ma-1) + this%emj(mj+m) * rcab(2,ma)
        cab(1,mj  ) =                    rcab(3,ma  )
        cab(2,mj  ) =                    rcab(4,ma  )
        
      else
        ma = ma+1
        mj = mj+1
        
        cab(1,mj) = this%emj(mj+m) * rcab(1,ma-1) + this%emj(mj+m+1) * rcab(1,ma)
        cab(2,mj) = this%emj(mj+m) * rcab(2,ma-1) + this%emj(mj+m+1) * rcab(2,ma)
      end if
    end do
    
    m = this%jmax
      !j == m
      ma = ma+1
      mj = mj+1
      
      cab(1,mj) = rcab(3,ma)
      cab(2,mj) = rcab(4,ma)
    
    !! Reindexing from order-fast jm to degree-fast mj indexing
    !! and synthethysing the real/imaginary parts into cmplx.
    m = 0
      !$omp simd
      do j = 0, this%jmax
        cjm(1+j*(j+1)/2)%re = cab(1,1+j)
        cjm(1+j*(j+1)/2)%im = 0._dbl
      end do
    
    do m = 1, this%jmax
      mj = m*(this%jmax+1)-m*(m+1)/2+1
      jm = m+1
      
      !$omp simd
      do j = m, this%jmax
        cjm(jm+j*(j+1)/2)%re = cab(1,mj+j)
        cjm(jm+j*(j+1)/2)%im = cab(2,mj+j)
      end do
    end do
    
    deallocate( cab )
    
  end procedure index_fwd_sub
  
end submodule r2c
