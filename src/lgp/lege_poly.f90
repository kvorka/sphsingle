module lege_poly
  use math
  implicit none
  
  type, public :: T_legep
    integer                             :: nFreq, jmax, jms, n, n_dbl, n_dbl_2, nrma
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
    
    procedure, public, pass :: is_rescale_sub
    procedure, public, pass :: bwd_indx_sub, index_bwd_sub, bwd_legesum_sub, bwd_sum1_sub, bwd_sum2_sub, bwd_shuffle_sub
    procedure, public, pass :: fwd_indx_sub, index_fwd_sub, fwd_legesum_sub, fwd_sum1_sub, fwd_sum2_sub, fwd_shuffle_sub
    
  end type T_legep
  
  real(kind=dbl), parameter :: deps  = 1.0d-15
  real(kind=qbl), parameter :: qeps  = 1.0d-28
  real(kind=qbl), parameter :: qpi   = acos(-1._qbl)
  real(kind=qbl), parameter :: qzero = 0._qbl
  
  interface
    module subroutine init_lege_sub(this, jmax, n, wfac)
      class(T_legep), intent(inout) :: this
      integer,        intent(in)    :: jmax, n
      real(kind=dbl), intent(in)    :: wfac
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
    
    module subroutine is_rescale_sub(this, rcab)
      class(T_legep), intent(in)    :: this
      real(kind=dbl), intent(inout) :: rcab(4,*)
    end subroutine is_rescale_sub
    
    module subroutine bwd_indx_sub(this, icab, ocab)
      class(T_legep), intent(in)  :: this
      real(kind=dbl), intent(in)  :: icab(2,*)
      real(kind=dbl), intent(out) :: ocab(2,2,*)
    end subroutine bwd_indx_sub
    
    module subroutine index_bwd_sub(this, cjm, rcab)
      class(T_legep),    intent(in)  :: this
      complex(kind=dbl), intent(in)  :: cjm(*)
      real(kind=dbl),    intent(out) :: rcab(*)
    end subroutine index_bwd_sub
    
    module subroutine bwd_sum1_sub(this, ma, pmm, pmj1, pmj, cc, swork)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: ma
      real(kind=dbl), intent(inout) :: pmj(ndbl,2,this%n_dbl_2), pmm(ndbl,2,this%n_dbl_2)
      real(kind=dbl), intent(in)    :: cc(4)
      real(kind=dbl), intent(out)   :: swork(ndbl,8,this%n_dbl_2), pmj1(ndbl,2,this%n_dbl_2)
    end subroutine bwd_sum1_sub
    
    module subroutine bwd_sum2_sub(this, ma, pmj1, pmj, cc, swork)
      class(T_legep), intent(in)     :: this
      integer,        intent(in)     :: ma
      real(kind=dbl), intent(inout)  :: pmj(ndbl,2,this%n_dbl_2)
      real(kind=dbl), intent(in)     :: cc(4), pmj1(ndbl,2,this%n_dbl_2)
      real(kind=dbl), intent(out)    :: swork(ndbl,8,this%n_dbl_2)
    end subroutine bwd_sum2_sub
    
    module subroutine bwd_shuffle_sub(this, swork, grid)
      class(T_legep), intent(in)  :: this
      real(kind=dbl), intent(in)  :: swork(ndbl,8,this%n_dbl_2)
      real(kind=dbl), intent(out) :: grid(ndbl,2,this%n_dbl_2,4)
    end subroutine bwd_shuffle_sub
    
    module subroutine bwd_legesum_sub(this, cc, grid, work)
      class(T_legep),         intent(in)  :: this
      real(kind=dbl),         intent(out) :: grid(4*this%n,0:this%jmax)
      real(kind=dbl), target, intent(out) :: work(*)
      real(kind=dbl),         intent(in)  :: cc(4,*)
    end subroutine bwd_legesum_sub
    
    module subroutine fwd_indx_sub(this, ocab, icab)
      class(T_legep), intent(in)  :: this
      real(kind=dbl), intent(in)  :: ocab(2,2,*)
      real(kind=dbl), intent(out) :: icab(2,*)
    end subroutine fwd_indx_sub
    
    module subroutine index_fwd_sub(this, rcab, cjm)
      class(T_legep),    intent(in)    :: this
      real(kind=dbl),    intent(inout) :: rcab(*)
      complex(kind=dbl), intent(out)   :: cjm(*)
    end subroutine index_fwd_sub
    
    module subroutine fwd_sum1_sub(this, ma, pmm, pmj1, pmj, swork, cr, acc, acc2)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: ma
      real(kind=dbl), intent(in)    :: swork(ndbl,8,this%n_dbl_2)
      real(kind=dbl), intent(inout) :: cr(4), pmj(ndbl,2,this%n_dbl_2), pmm(ndbl,2,this%n_dbl_2)
      real(kind=dbl), intent(out)   :: acc(ndbl,4), acc2(ndbl,4), pmj1(ndbl,2,this%n_dbl_2)
    end subroutine fwd_sum1_sub
    
    module subroutine fwd_sum2_sub(this, ima, pmj1, pmj, swork, cr, acc, acc2)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: ima
      real(kind=dbl), intent(inout) :: pmj(ndbl,2,this%n_dbl_2)
      real(kind=dbl), intent(in)    :: swork(ndbl,8,this%n_dbl_2), pmj1(ndbl,2,this%n_dbl_2)
      real(kind=dbl), intent(inout) :: cr(4)
      real(kind=dbl), intent(out)   :: acc(ndbl,4), acc2(ndbl,4)
    end subroutine fwd_sum2_sub
    
    module subroutine fwd_shuffle_sub(this, grid, swork)
      class(T_legep), intent(in)  :: this
      real(kind=dbl), intent(in)  :: grid(ndbl,2,this%n_dbl_2,4)
      real(kind=dbl), intent(out) :: swork(ndbl,8,this%n_dbl_2)
    end subroutine fwd_shuffle_sub
    
    module subroutine fwd_legesum_sub(this, grid, cr, work)
      class(T_legep),         intent(in)    :: this
      real(kind=dbl),         intent(inout) :: grid(4*this%n,0:this%jmax)
      real(kind=dbl),         intent(inout) :: cr(4,*)
      real(kind=dbl), target, intent(out)   :: work(*)
    end subroutine fwd_legesum_sub
  end interface
  
end module lege_poly
