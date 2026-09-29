submodule (lege_poly) c2r
  implicit none; contains
  
  module procedure index_bwd_sub
    integer                     :: m, j, jm, mj, ma
    real(kind=dbl), allocatable :: cab(:,:)
    
    !! Reindexing of real/imaginary parts of the original sequence 
    !! from order-fast jm to degree-fast mj indexing.
    allocate( cab(2,this%jms) )
    
    do m = 0, this%jmax
      mj = m*(this%jmax+1)-m*(m+1)/2+1
      
      !$omp simd
      do j = m, this%jmax
        jm = j*(j+1)/2+m+1
        
        cab(1,mj+j) = cjm(jm)%re
        cab(2,mj+j) = cjm(jm)%im
      end do
    end do
    
    !! Reindexing for transform. Quadruplets of rescaled coefficients
    !! are prepared. These include real/imaginary, odd degree/even degree
    !! components, respectively, for easier caching.
    m = 0
      !j == m
        ma = 1
        mj = 1
        
        rcab(1,ma) = cab(1,mj+1) * this%emj(mj+1)
        rcab(2,ma) = cab(2,mj+1) * this%emj(mj+1)
        rcab(3,ma) = cab(1,mj)
        rcab(4,ma) = cab(2,mj)
      
      !$omp simd
      do j = 1, (this%jmax-1)/2
        ma = ma+1
        mj = mj+2
        
        rcab(1,ma) = this%emj(mj) * cab(1,mj-1) + this%emj(mj+1) * cab(1,mj+1)
        rcab(2,ma) = this%emj(mj) * cab(2,mj-1) + this%emj(mj+1) * cab(2,mj+1)
        rcab(3,ma) =                cab(1,mj)
        rcab(4,ma) =                cab(2,mj)
      end do
      
      !j == jmax
      if ( mod(this%jmax,2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        rcab(1,ma) = this%emj(mj) * cab(1,mj-1)
        rcab(2,ma) = this%emj(mj) * cab(2,mj-1)
        rcab(3,ma) =                cab(1,mj)
        rcab(4,ma) =                cab(2,mj)
        
      else
        ma = ma+1
        mj = mj+1
        
        rcab(1,ma) = this%emj(mj+1) * cab(1,mj)
        rcab(2,ma) = this%emj(mj+1) * cab(2,mj)
        rcab(3,ma) = 0._dbl
        rcab(4,ma) = 0._dbl
      end if
    
    do m = 1, this%jmax-1
      !j == m
        ma = ma+1
        mj = mj+1
        
        rcab(1,ma) = cab(1,mj+1) * this%emj(mj+m+1)
        rcab(2,ma) = cab(2,mj+1) * this%emj(mj+m+1)
        rcab(3,ma) = cab(1,mj)
        rcab(4,ma) = cab(2,mj)
      
      !$omp simd
      do j = 1, (this%jmax-1-m)/2
        ma = ma+1
        mj = mj+2
        
        rcab(1,ma) = this%emj(mj+m) * cab(1,mj-1) + this%emj(mj+m+1) * cab(1,mj+1)
        rcab(2,ma) = this%emj(mj+m) * cab(2,mj-1) + this%emj(mj+m+1) * cab(2,mj+1)
        rcab(3,ma) =                  cab(1,mj)
        rcab(4,ma) =                  cab(2,mj)
      end do
      
      !j == jmax
      if ( mod((this%jmax-m),2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        rcab(1,ma) = this%emj(mj+m) * cab(1,mj-1)
        rcab(2,ma) = this%emj(mj+m) * cab(2,mj-1)
        rcab(3,ma) =                  cab(1,mj)
        rcab(4,ma) =                  cab(2,mj)
      
      else
        ma = ma+1
        mj = mj+1
        
        rcab(1,ma) = this%emj(mj+m+1) * cab(1,mj)
        rcab(2,ma) = this%emj(mj+m+1) * cab(2,mj)
        rcab(3,ma) = 0._dbl
        rcab(4,ma) = 0._dbl
      end if
    end do
    
    m = this%jmax
      !j == m
        ma = ma+1
        mj = mj+1
        
        rcab(1,ma) = 0._dbl
        rcab(2,ma) = 0._dbl
        rcab(3,ma) = cab(1,mj)
        rcab(4,ma) = cab(2,mj)
    
    deallocate( cab )
    
    !! Last rescale required by the Ishioka recursion. This call is not
    !! inlined because the same rescale is required after the forward
    !! transform. Vectorized inside.
    call is_rescale_sub( this%nrma, this%amj, rcab )
    
  end procedure index_bwd_sub
  
end submodule c2r