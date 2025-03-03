submodule (lege_poly) c2r
  implicit none; contains
  
  module procedure index_bwd_sub
    integer                     :: m, j, jm, mj, ma
    real(kind=dbl), allocatable :: cab(:,:)
    
    allocate( cab(2,this%jms) )
    
    do m = 0, this%jmax
      do j = m, this%jmax
        mj = m*(this%jmax+1)-m*(m+1)/2+j+1
        jm = j*(j+1)/2+m+1
        
        cab(1,mj) = cjm(jm)%re
        cab(2,mj) = cjm(jm)%im
      end do
    end do
    
    m = 0
      !j == m
        ma = 1
        mj = 1
        
        rcab(1,1,ma) = cab(1,mj+1) * this%emj(mj+1)
        rcab(2,1,ma) = cab(2,mj+1) * this%emj(mj+1)
        rcab(1,2,ma) = cab(1,mj)
        rcab(2,2,ma) = cab(2,mj)
      
      do j = 1, (this%jmax-1)/2
        ma = ma+1
        mj = mj+2
        
        rcab(1,1,ma) = this%emj(mj) * cab(1,mj-1) + this%emj(mj+1) * cab(1,mj+1)
        rcab(2,1,ma) = this%emj(mj) * cab(2,mj-1) + this%emj(mj+1) * cab(2,mj+1)
        rcab(1,2,ma) =                cab(1,mj)
        rcab(2,2,ma) =                cab(2,mj)
      end do
      
      !j == this%jmax
      if ( mod((this%jmax),2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        rcab(1,1,ma) = this%emj(mj) * cab(1,mj-1)
        rcab(2,1,ma) = this%emj(mj) * cab(2,mj-1)
        rcab(1,2,ma) =                cab(1,mj)
        rcab(2,2,ma) =                cab(2,mj)
        
      else
        ma = ma+1
        mj = mj+1
        
        rcab(1,1,ma) = this%emj(mj+1) * cab(1,mj)
        rcab(2,1,ma) = this%emj(mj+1) * cab(2,mj)
      end if
    
    do m = 1, this%jmax-1
      !j == m
        ma = ma+1
        mj = mj+1
        
        rcab(1,1,ma) = cab(1,mj+1) * this%emj(mj+m+1)
        rcab(2,1,ma) = cab(2,mj+1) * this%emj(mj+m+1)
        rcab(1,2,ma) = cab(1,mj)
        rcab(2,2,ma) = cab(2,mj)
      
      do j = 1, (this%jmax-1-m)/2
        ma = ma+1
        mj = mj+2
        
        rcab(1,1,ma) = this%emj(mj+m) * cab(1,mj-1) + this%emj(mj+m+1) * cab(1,mj+1)
        rcab(2,1,ma) = this%emj(mj+m) * cab(2,mj-1) + this%emj(mj+m+1) * cab(2,mj+1)
        rcab(1,2,ma) =                  cab(1,mj)
        rcab(2,2,ma) =                  cab(2,mj)
      end do
      
      !j == this%jmax
      if ( mod((this%jmax-m),2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        rcab(1,1,ma) = this%emj(mj+m) * cab(1,mj-1)
        rcab(2,1,ma) = this%emj(mj+m) * cab(2,mj-1)
        rcab(1,2,ma) =                  cab(1,mj)
        rcab(2,2,ma) =                  cab(2,mj)
      
      else
        ma = ma+1
        mj = mj+1
        
        rcab(1,1,ma) = this%emj(mj+m+1) * cab(1,mj)
        rcab(2,1,ma) = this%emj(mj+m+1) * cab(2,mj)
      end if
    end do
    
    m = this%jmax
      !j == m
        ma = ma+1
        mj = mj+1
        
        rcab(1,2,ma) = cab(1,mj)
        rcab(2,2,ma) = cab(2,mj)
    
    deallocate( cab )
    
    call this%is_rescale_sub( rcab )
    
  end procedure index_bwd_sub
  
end submodule c2r