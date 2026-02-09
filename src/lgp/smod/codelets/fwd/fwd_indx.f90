submodule(lege_poly) fwd_indx
  implicit none; contains
  
  module procedure fwd_indx_sub
    integer :: m, j, mj, ma
    
    m = 0
      !j == m
        ma = 1
        mj = 1
        
        icab(1,mj) = ocab(1,2,ma)
        icab(2,mj) = ocab(2,2,ma)
      
      do j = 1, (this%jmax-1)/2
        ma = ma+1
        mj = mj+2
        
        icab(1,mj-1) = this%emj(mj-1) * ocab(1,1,ma-1) + this%emj(mj) * ocab(1,1,ma)
        icab(2,mj-1) = this%emj(mj-1) * ocab(2,1,ma-1) + this%emj(mj) * ocab(2,1,ma)
        icab(1,mj  ) =                  ocab(1,2,ma  )
        icab(2,mj  ) =                  ocab(2,2,ma  )
      end do
      
      !j == jmax
      if ( mod((this%jmax),2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        icab(1,mj-1) = this%emj(mj-1) * ocab(1,1,ma-1) + this%emj(mj) * ocab(1,1,ma)
        icab(2,mj-1) = this%emj(mj-1) * ocab(2,1,ma-1) + this%emj(mj) * ocab(2,1,ma)
        icab(1,mj  ) =                  ocab(1,2,ma  )
        icab(2,mj  ) =                  ocab(2,2,ma  )
      
      else
        ma = ma+1
        mj = mj+1
        
        icab(1,mj) = this%emj(mj) * ocab(1,1,ma-1) + this%emj(mj+1) * ocab(1,1,ma)
        icab(2,mj) = this%emj(mj) * ocab(2,1,ma-1) + this%emj(mj+1) * ocab(2,1,ma)
      end if
    
    do m = 1, this%jmax-1
      !j == m
        ma = ma+1
        mj = mj+1
        
        icab(1,mj) = ocab(1,2,ma)
        icab(2,mj) = ocab(2,2,ma)
      
      do j = 1, (this%jmax-m-1)/2
        ma = ma+1
        mj = mj+2
        
        icab(1,mj-1) = this%emj(mj+m-1) * ocab(1,1,ma-1) + this%emj(mj+m) * ocab(1,1,ma)
        icab(2,mj-1) = this%emj(mj+m-1) * ocab(2,1,ma-1) + this%emj(mj+m) * ocab(2,1,ma)
        icab(1,mj  ) =                    ocab(1,2,ma  )
        icab(2,mj  ) =                    ocab(2,2,ma  )
      end do
      
      !j == jmax
      if ( mod((this%jmax-m),2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        icab(1,mj-1) = this%emj(mj+m-1) * ocab(1,1,ma-1) + this%emj(mj+m) * ocab(1,1,ma)
        icab(2,mj-1) = this%emj(mj+m-1) * ocab(2,1,ma-1) + this%emj(mj+m) * ocab(2,1,ma)
        icab(1,mj  ) =                    ocab(1,2,ma  )
        icab(2,mj  ) =                    ocab(2,2,ma  )
        
      else
        ma = ma+1
        mj = mj+1
        
        icab(1,mj) = this%emj(mj+m) * ocab(1,1,ma-1) + this%emj(mj+m+1) * ocab(1,1,ma)
        icab(2,mj) = this%emj(mj+m) * ocab(2,1,ma-1) + this%emj(mj+m+1) * ocab(2,1,ma)
      end if
    end do
    
    m = this%jmax
      !j == m
      ma = ma+1
      mj = mj+1
      
      icab(1,mj) = ocab(1,2,ma)
      icab(2,mj) = ocab(2,2,ma)
  
  end procedure fwd_indx_sub
  
end submodule fwd_indx