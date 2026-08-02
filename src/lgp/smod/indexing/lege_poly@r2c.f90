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
        
        cab(1,mj) = rcab(1,2,ma)
        cab(2,mj) = rcab(2,2,ma)
      
      !$omp simd
      do j = 1, (this%jmax-1)/2
        ma = ma+1
        mj = mj+2
        
        cab(1,mj-1) = this%emj(mj-1) * rcab(1,1,ma-1) + this%emj(mj) * rcab(1,1,ma)
        cab(2,mj-1) = this%emj(mj-1) * rcab(2,1,ma-1) + this%emj(mj) * rcab(2,1,ma)
        cab(1,mj  ) =                  rcab(1,2,ma  )
        cab(2,mj  ) =                  rcab(2,2,ma  )
      end do
      
      !j == jmax
      if ( mod(this%jmax,2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        cab(1,mj-1) = this%emj(mj-1) * rcab(1,1,ma-1) + this%emj(mj) * rcab(1,1,ma)
        cab(2,mj-1) = this%emj(mj-1) * rcab(2,1,ma-1) + this%emj(mj) * rcab(2,1,ma)
        cab(1,mj  ) =                  rcab(1,2,ma  )
        cab(2,mj  ) =                  rcab(2,2,ma  )
      
      else
        ma = ma+1
        mj = mj+1
        
        cab(1,mj) = this%emj(mj) * rcab(1,1,ma-1) + this%emj(mj+1) * rcab(1,1,ma)
        cab(2,mj) = this%emj(mj) * rcab(2,1,ma-1) + this%emj(mj+1) * rcab(2,1,ma)
      end if
    
    do m = 1, this%jmax-1
      !j == m
        ma = ma+1
        mj = mj+1
        
        cab(1,mj) = rcab(1,2,ma)
        cab(2,mj) = rcab(2,2,ma)
      
      !$omp simd
      do j = 1, (this%jmax-m-1)/2
        ma = ma+1
        mj = mj+2
        
        cab(1,mj-1) = this%emj(mj+m-1) * rcab(1,1,ma-1) + this%emj(mj+m) * rcab(1,1,ma)
        cab(2,mj-1) = this%emj(mj+m-1) * rcab(2,1,ma-1) + this%emj(mj+m) * rcab(2,1,ma)
        cab(1,mj  ) =                    rcab(1,2,ma  )
        cab(2,mj  ) =                    rcab(2,2,ma  )
      end do
      
      !j == jmax
      if ( mod((this%jmax-m),2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        cab(1,mj-1) = this%emj(mj+m-1) * rcab(1,1,ma-1) + this%emj(mj+m) * rcab(1,1,ma)
        cab(2,mj-1) = this%emj(mj+m-1) * rcab(2,1,ma-1) + this%emj(mj+m) * rcab(2,1,ma)
        cab(1,mj  ) =                    rcab(1,2,ma  )
        cab(2,mj  ) =                    rcab(2,2,ma  )
        
      else
        ma = ma+1
        mj = mj+1
        
        cab(1,mj) = this%emj(mj+m) * rcab(1,1,ma-1) + this%emj(mj+m+1) * rcab(1,1,ma)
        cab(2,mj) = this%emj(mj+m) * rcab(2,1,ma-1) + this%emj(mj+m+1) * rcab(2,1,ma)
      end if
    end do
    
    m = this%jmax
      !j == m
      ma = ma+1
      mj = mj+1
      
      cab(1,mj) = rcab(1,2,ma)
      cab(2,mj) = rcab(2,2,ma)
    
    !! Reindexing from order-fast jm to degree-fast mj indexing
    !! and synthethysing the real/imaginary parts into cmplx.
    m = 0
      !$omp simd
      do j = 0, this%jmax
        jm = j*(j+1)/2+1
        
        cjm(jm)%re = cab(1,1+j)
        cjm(jm)%im = zero
      end do
    
    !GCC$ unroll 4
    !DIR$ unroll (4)
    do m = 1, this%jmax
      mj = m*(this%jmax+1)-m*(m+1)/2+1
      
      !$omp simd
      do j = m, this%jmax
        jm = j*(j+1)/2+m+1
        
        cjm(jm)%re = cab(1,mj+j)
        cjm(jm)%im = cab(2,mj+j)
      end do
    end do
    
    deallocate( cab )
    
  end procedure index_fwd_sub
  
end submodule r2c
