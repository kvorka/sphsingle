module lege_poly
  use math
  implicit none
  
  type, public :: T_legep
    integer                             :: nFreq, jmax, jms, n, n_dbl, nrma
    integer,        allocatable         :: mamj(:)
    real(kind=dbl), allocatable         :: emj(:), fmj(:,:), amj(:)
    real(kind=dbl), pointer, contiguous :: cosx(:), sinx(:), cosx2(:), wght(:)
    type(c_ptr)                         :: c_cosx, c_sinx, c_cosx2, c_wght
    
    contains
    
    procedure, public,  pass :: init_sub       => init_lege_sub
    procedure, private, pass :: roots_sub      => find_roots_sub
    procedure, private, pass :: coeffs_sub     => compute_coeffs_sub
    procedure, private, pass :: get_nma_sub    => get_nma_sub
    procedure, public,  pass :: deallocate_sub => deallocate_lege_sub
    
    procedure, public, pass :: index_bwd_sub, bwd_legesum_sub
    procedure, public, pass :: index_fwd_sub, fwd_legesum_sub
    
  end type T_legep
  
  real(kind=dbl), parameter :: deps  = 1.0d-15
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
    
    module subroutine get_nma_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine get_nma_sub
    
    module subroutine index_bwd_sub(this, cjm, rcab)
      class(T_legep),    intent(in)  :: this
      complex(kind=dbl), intent(in)  :: cjm(*)
      real(kind=dbl),    intent(out) :: rcab(*)
    end subroutine index_bwd_sub
    
    module subroutine bwd_legesum_sub(this, cc, grid, work)
      class(T_legep),         intent(in)  :: this
      real(kind=dbl),         intent(out) :: grid(this%n,4,0:this%jmax)
      real(kind=dbl), target, intent(out) :: work(*)
      real(kind=dbl),         intent(in)  :: cc(4,*)
    end subroutine bwd_legesum_sub
    
    module subroutine index_fwd_sub(this, rcab, cjm)
      class(T_legep),    intent(in)    :: this
      real(kind=dbl),    intent(inout) :: rcab(*)
      complex(kind=dbl), intent(out)   :: cjm(*)
    end subroutine index_fwd_sub
    
    module subroutine fwd_legesum_sub(this, grid, cr, work)
      class(T_legep),         intent(in)    :: this
      real(kind=dbl),         intent(inout) :: grid(this%n,4,0:this%jmax)
      real(kind=dbl),         intent(inout) :: cr(4,*)
      real(kind=dbl), target, intent(out)   :: work(*)
    end subroutine fwd_legesum_sub
  end interface
  
  interface
    module subroutine is_rescale_sub(nrma, amj, rcab)
      integer,        intent(in)    :: nrma
      real(kind=dbl), intent(in)    :: amj(nrma)
      real(kind=dbl), intent(inout) :: rcab(4,nrma)
    end subroutine is_rescale_sub
    
    module subroutine bwd_indx_sub(jmax, emj, icab, ocab)
      integer,        intent(in)  :: jmax
      real(kind=dbl), intent(in)  :: icab(2,*), emj(*)
      real(kind=dbl), intent(out) :: ocab(2,2,*)
    end subroutine bwd_indx_sub
    
    module subroutine bwd_sum1_sub(n1, ma, fmj, sinx, cosx, pmm, pmj1, pmj, cc, s1, s2, s3, s4)
      integer,         intent(in)    :: n1, ma
      real(kind=dbl),  intent(inout) :: pmj(ndbl,n1), pmm(ndbl,n1)
      real(kind=dbl),  intent(in)    :: cc(4), fmj, sinx(ndbl,n1), cosx(ndbl,n1)
      real(kind=dbl),  intent(out)   :: pmj1(ndbl,n1), s1(ndbl,n1), s2(ndbl,n1), s3(ndbl,n1), s4(ndbl,n1)
    end subroutine bwd_sum1_sub
    
    module subroutine bwd_sum2_sub(n1, fmj, cosx2, pmj1, pmj, cc, s1, s2, s3, s4)
      integer,        intent(in)    :: n1
      real(kind=dbl), intent(inout) :: pmj(ndbl,n1)
      real(kind=dbl), intent(in)    :: fmj(2), cosx2(ndbl,n1), cc(4), pmj1(ndbl,n1)
      real(kind=dbl), intent(out)   :: s1(ndbl,n1), s2(ndbl,n1), s3(ndbl,n1), s4(ndbl,n1)
    end subroutine bwd_sum2_sub
    
    module subroutine bwd_shuffle_sub(n1, cosx, g1, g2, g3, g4)
      integer,        intent(in)    :: n1
      real(kind=dbl), intent(in)    :: cosx(ndbl,n1)
      real(kind=dbl), intent(inout) :: g1(ndbl,n1), g2(ndbl,n1), g3(ndbl,n1), g4(ndbl,n1)
    end subroutine bwd_shuffle_sub
    
    module subroutine fwd_shuffle_sub(n1, cosx, wght, g1, g2, g3, g4)
      integer,        intent(in)    :: n1
      real(kind=dbl), intent(in)    :: cosx(ndbl,n1), wght(ndbl,n1)
      real(kind=dbl), intent(inout) :: g1(ndbl,n1), g2(ndbl,n1), g3(ndbl,n1), g4(ndbl,n1)
    end subroutine fwd_shuffle_sub
    
    module subroutine fwd_sum1_sub(n1, ma, fmj, sinx, cosx, pmm, pmj1, pmj, s1, s2, s3, s4, cr, acc, acc2)
      integer,        intent(in)    :: n1, ma
      real(kind=dbl), intent(in)    :: s1(ndbl,n1), s2(ndbl,n1), s3(ndbl,n1), s4(ndbl,n1), fmj, sinx(ndbl,n1), cosx(ndbl,n1)
      real(kind=dbl), intent(inout) :: cr(4), pmj(ndbl,n1), pmm(ndbl,n1)
      real(kind=dbl), intent(out)   :: acc(ndbl,4), acc2(ndbl,4), pmj1(ndbl,n1)
    end subroutine fwd_sum1_sub
    
    module subroutine fwd_sum2_sub(n1, fmj, cosx2, pmj1, pmj, s1, s2, s3, s4, cr, acc, acc2)
      integer,        intent(in)    :: n1
      real(kind=dbl), intent(inout) :: pmj(ndbl,n1)
      real(kind=dbl), intent(in)    :: s1(ndbl,n1), s2(ndbl,n1), s3(ndbl,n1), s4(ndbl,n1), pmj1(ndbl,n1), fmj(2), cosx2(ndbl,n1)
      real(kind=dbl), intent(inout) :: cr(4)
      real(kind=dbl), intent(out)   :: acc(ndbl,4), acc2(ndbl,4)
    end subroutine fwd_sum2_sub
    
    module subroutine fwd_indx_sub(jmax, emj, ocab, icab)
      integer,        intent(in)  :: jmax
      real(kind=dbl), intent(in)  :: ocab(2,2,*), emj(*)
      real(kind=dbl), intent(out) :: icab(2,*)
    end subroutine fwd_indx_sub
  end interface
  
end module lege_poly
