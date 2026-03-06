submodule (lege_poly) get_nma
  implicit none; contains
  
  module procedure get_nma_sub
    integer :: j, m, ma
    
    this%nrma = 0
      do m = 0, this%jmax
        this%nrma = this%nrma+1
        
        if ( m < this%jmax) then
          do j = 1, (this%jmax-1-m)/2
            this%nrma = this%nrma+1
          end do
          
          this%nrma = this%nrma+1
        end if
      end do
    
    allocate( this%mamj(0:this%jmax+1) )
    
    ma = 0
        
    do m = 0, this%jmax
      !j = m
        ma = ma+1
        this%mamj(m) = ma
      
      do j = 1, (this%jmax-m)/2
        ma = ma+1
      end do
      
      if ( mod((this%jmax-m),2) /= 0 ) then
        ma = ma+1
      end if
    end do
    
    this%mamj(this%jmax+1) = ma+1
    
  end procedure get_nma_sub

end submodule get_nma