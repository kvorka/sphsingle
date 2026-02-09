submodule(lege_poly) bwd_indx
  implicit none; contains
  
  module procedure bwd_indx_sub
    integer :: m, j, mj, ma
    
    m = 0
      !j == m
        ma = 1
        mj = 1
        
        ocab(1,1,ma) = icab(1,mj+1) * this%emj(mj+1)
        ocab(2,1,ma) = icab(2,mj+1) * this%emj(mj+1)
        ocab(1,2,ma) = icab(1,mj)
        ocab(2,2,ma) = icab(2,mj)
      
      do j = 1, (this%jmax-1)/2
        ma = ma+1
        mj = mj+2
        
        ocab(1,1,ma) = this%emj(mj) * icab(1,mj-1) + this%emj(mj+1) * icab(1,mj+1)
        ocab(2,1,ma) = this%emj(mj) * icab(2,mj-1) + this%emj(mj+1) * icab(2,mj+1)
        ocab(1,2,ma) =                icab(1,mj)
        ocab(2,2,ma) =                icab(2,mj)
      end do
      
      !j == jmax
      if ( mod(this%jmax,2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        ocab(1,1,ma) = this%emj(mj) * icab(1,mj-1)
        ocab(2,1,ma) = this%emj(mj) * icab(2,mj-1)
        ocab(1,2,ma) =                icab(1,mj)
        ocab(2,2,ma) =                icab(2,mj)
        
      else
        ma = ma+1
        mj = mj+1
        
        ocab(1,1,ma) = this%emj(mj+1) * icab(1,mj)
        ocab(2,1,ma) = this%emj(mj+1) * icab(2,mj)
      end if
    
    do m = 1, this%jmax-1
      !j == m
        ma = ma+1
        mj = mj+1
        
        ocab(1,1,ma) = icab(1,mj+1) * this%emj(mj+m+1)
        ocab(2,1,ma) = icab(2,mj+1) * this%emj(mj+m+1)
        ocab(1,2,ma) = icab(1,mj)
        ocab(2,2,ma) = icab(2,mj)
      
      do j = 1, (this%jmax-1-m)/2
        ma = ma+1
        mj = mj+2
        
        ocab(1,1,ma) = this%emj(mj+m) * icab(1,mj-1) + this%emj(mj+m+1) * icab(1,mj+1)
        ocab(2,1,ma) = this%emj(mj+m) * icab(2,mj-1) + this%emj(mj+m+1) * icab(2,mj+1)
        ocab(1,2,ma) =                  icab(1,mj)
        ocab(2,2,ma) =                  icab(2,mj)
      end do
      
      !j == jmax
      if ( mod((this%jmax-m),2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        ocab(1,1,ma) = this%emj(mj+m) * icab(1,mj-1)
        ocab(2,1,ma) = this%emj(mj+m) * icab(2,mj-1)
        ocab(1,2,ma) =                  icab(1,mj)
        ocab(2,2,ma) =                  icab(2,mj)
      
      else
        ma = ma+1
        mj = mj+1
        
        ocab(1,1,ma) = this%emj(mj+m+1) * icab(1,mj)
        ocab(2,1,ma) = this%emj(mj+m+1) * icab(2,mj)
      end if
    end do
    
    m = this%jmax
      !j == m
        ma = ma+1
        mj = mj+1
        
        ocab(1,2,ma) = icab(1,mj)
        ocab(2,2,ma) = icab(2,mj)
  
  end procedure bwd_indx_sub
  
end submodule bwd_indx