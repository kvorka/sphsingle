module lege_poly
  use math
  implicit none
  
  type, public :: T_legep
    integer                             :: nFreq, jmax, jms, n, n_dbl, nrma
    integer,        allocatable         :: mamj(:)
    real(kind=dbl), allocatable         :: emj(:), fmj(:,:), amj(:)
    real(kind=dbl), pointer, contiguous :: cosx(:), cosx2(:), wght(:), pmm(:,:,:)
    type(c_ptr)                         :: c_cosx, c_cosx2, c_wght, c_pmm
    
    contains
    
    procedure, public,  pass :: init_sub       => init_lege_sub
    procedure, public,  pass :: deallocate_sub => deallocate_lege_sub
    
    procedure, public, pass :: index_bwd_sub, bwd_legesum_sub
    procedure, public, pass :: index_fwd_sub, fwd_legesum_sub
    
  end type T_legep
  
  real(kind=qbl), parameter :: qeps  = 1.0d-28
  real(kind=qbl), parameter :: qpi   = acos(-1._qbl)
  real(kind=qbl), parameter :: qzero = 0._qbl
  
  interface
    module subroutine init_lege_sub(this, jmax, wfac)
      class(T_legep), intent(inout) :: this
      integer,        intent(in)    :: jmax
      integer,        intent(in)    :: wfac
    end subroutine init_lege_sub
    
    module subroutine deallocate_lege_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine deallocate_lege_sub
    
    module subroutine find_roots_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine find_roots_sub
    
    module subroutine compute_coeffs_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine compute_coeffs_sub
    
    module subroutine compute_pmm_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine compute_pmm_sub
    
    module subroutine index_bwd_sub(this, cjm, rcab)
      class(T_legep),    intent(in)  :: this
      complex(kind=dbl), intent(in)  :: cjm(*)
      real(kind=dbl),    intent(out) :: rcab(2,2,*)
    end subroutine index_bwd_sub
    
    module subroutine bwd_legesum_sub(this, cc, grid)
      class(T_legep), intent(in)  :: this
      real(kind=dbl), intent(out) :: grid(4*this%n,0:this%jmax)
      real(kind=dbl), intent(in)  :: cc(4,*)
    end subroutine bwd_legesum_sub
    
    module subroutine index_fwd_sub(this, rcab, cjm)
      class(T_legep),    intent(in)    :: this
      real(kind=dbl),    intent(inout) :: rcab(2,2,*)
      complex(kind=dbl), intent(out)   :: cjm(*)
    end subroutine index_fwd_sub
    
    module subroutine fwd_legesum_sub(this, grid, cr)
      class(T_legep), intent(in)    :: this
      real(kind=dbl), intent(inout) :: grid(4*this%n,0:this%jmax)
      real(kind=dbl), intent(inout) :: cr(4,*)
    end subroutine fwd_legesum_sub
  end interface
  
  interface
    module subroutine is_rescale_sub(nrma, amj, rcab)
      integer,        intent(in)    :: nrma
      real(kind=dbl), intent(in)    :: amj(nrma)
      real(kind=dbl), intent(inout) :: rcab(4,nrma)
    end subroutine is_rescale_sub
    
    module subroutine bwd_sum_m_sub(n1, ma1, ma2, fmj, cosx, cosx2, pmm, pmj1, pmj, cc, swork, grid)
      integer,                             intent(in)    :: n1, ma1, ma2
      real(kind=dbl),                      intent(in)    :: fmj(2,ma1:ma2), cosx(ndbl,n1), cosx2(ndbl,n1)
      real(kind=dbl),                      intent(in)    :: cc(4,ma1:ma2)
      real(kind=dbl),                      intent(out)   :: swork(ndbl,4,n1), grid(ndbl,n1,4)
      real(kind=dbl), pointer, contiguous, intent(in)    :: pmm(:,:)
      real(kind=dbl), pointer, contiguous, intent(inout) :: pmj(:,:), pmj1(:,:)
    end subroutine bwd_sum_m_sub
    
    module subroutine fwd_sum_m_sub(n1, ma1, ma2, fmj, cosx, cosx2, wght, pmm, pmj1, pmj, swork, cr, acc, acc2, grid)
      integer,        intent(in)                         :: n1, ma1, ma2
      real(kind=dbl), intent(in)                         :: fmj(2,ma1:ma2), cosx(ndbl,n1), cosx2(ndbl,n1), &
                                                          & wght(ndbl,n1), grid(ndbl,n1,4)
      real(kind=dbl), intent(inout)                      :: cr(4,ma1:ma2)
      real(kind=dbl), intent(out)                        :: acc(ndbl,4), acc2(ndbl,4), swork(ndbl,4,n1)
      real(kind=dbl), pointer, contiguous, intent(in)    :: pmm(:,:)
      real(kind=dbl), pointer, contiguous, intent(inout) :: pmj(:,:), pmj1(:,:)
    end subroutine fwd_sum_m_sub
  end interface
  
end module lege_poly
