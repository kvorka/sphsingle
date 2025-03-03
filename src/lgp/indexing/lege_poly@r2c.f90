submodule (lege_poly) r2c
  implicit none; contains
  
  module procedure index_fwd_sub
    integer                     :: m, j, jm, mj, ma
    real(kind=dbl), allocatable :: cab(:,:)
    
    call this%is_rescale_sub( rcab )
    
    allocate( cab(2,this%jms) )
    
    m = 0
      !j == m
        ma = 1
        mj = 1
        
        cab(1,mj) = rcab(1,2,ma)
        cab(2,mj) = rcab(2,2,ma)
      
      do j = 1, (this%jmax-1)/2
        ma = ma+1
        mj = mj+2
        
        cab(1,mj-1) = this%emj(mj-1) * rcab(1,1,ma-1) + this%emj(mj) * rcab(1,1,ma)
        cab(2,mj-1) = this%emj(mj-1) * rcab(2,1,ma-1) + this%emj(mj) * rcab(2,1,ma)
        cab(1,mj  ) =                  rcab(1,2,ma  )
        cab(2,mj  ) =                  rcab(2,2,ma  )
      end do
      
      !j == this%jmax
      if ( mod((this%jmax),2) == 0 ) then
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
      
      do j = 1, (this%jmax-m-1)/2
        ma = ma+1
        mj = mj+2
        
        cab(1,mj-1) = this%emj(mj+m-1) * rcab(1,1,ma-1) + this%emj(mj+m) * rcab(1,1,ma)
        cab(2,mj-1) = this%emj(mj+m-1) * rcab(2,1,ma-1) + this%emj(mj+m) * rcab(2,1,ma)
        cab(1,mj  ) =                    rcab(1,2,ma  )
        cab(2,mj  ) =                    rcab(2,2,ma  )
      end do
      
      !j == this%jmax
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
    
    do j = 0, this%jmax
      m = 0
        mj = m*(this%jmax+1)-m*(m+1)/2+j+1
        jm = j*(j+1)/2+m+1
        
        cjm(jm) = cmplx( cab(1,mj), 0._dbl, kind=dbl )
      
      do m = 1, j
        mj = m*(this%jmax+1)-m*(m+1)/2+j+1
        jm = j*(j+1)/2+m+1
        
        cjm(jm) = cmplx( cab(1,mj), cab(2,mj), kind=dbl )
      end do
    end do
    
    deallocate( cab )
    
  end procedure index_fwd_sub
  
end submodule r2c
