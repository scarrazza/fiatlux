!----------------------------------------------------------------------
! A module to define all parameters of the run, which can either be read from
! an input card (powheg.input, vbfnlo.input) or from command line arguments
module parameters
  use sub_defs_io; use types; use consts_dp
  !use integration
  use types
  implicit none

  private
  !real(dp), parameter, public :: gfermi = 1.16637e-5_dp
  !real(dp), parameter, public :: gev2pb = 389379660.0_dp
  !real(dp), parameter, public :: eps    = 1.0e-14_dp
  !real(dp), public :: mh, mh_sq, hwidth, mt, mb
  !real(dp), public :: sqrts, S, Q0_cut_sq
  !real(dp), public :: hmasswindow, toyas
  !integer,  public :: order_min, order_max!, scale_choice
  !logical,  public :: higgs_use_BW, toypdf
  !character * 4, public :: seedstr    

  ! for now we fix the number of flavours -- to be
  ! revisited later...
  integer,  public :: nflav = 5
  
  real(dp), public :: xmuf=one, xmur = one, Qmin = one

  ! if test_Q0 > zero, then PDF is read at that scale
  ! and then evolved up with muR/Q = muR_pdf
  real(dp), public :: test_Q0 = -one , muR_PDF = one

  ! some electroweak parameters
  real(dp), public, parameter :: mw = 80.398_dp, mz = 91.187_dp
  real(dp), public, parameter :: w_width = 2.141_dp, z_width = 2.4952_dp
  real(dp), public, parameter :: sin_thw = 1.0_dp - (mw/mz)**2

  !public :: set_parameters

contains

!   ! set all parameters from input card or command line arguments
!   subroutine set_parameters (read_from_card)
!     logical, intent(in) :: read_from_card
!     real(dp) :: powheginput, vbfnloinput
!     real(dp) :: ebeam1, ebeam2
!     integer  :: idum, iseed
!     common/ranno/idum
!     character * 6 WHCPRG
!     common/cWHCPRG/WHCPRG
! 
!     if (read_from_card) then
!        ! if read_from_card is true, read all parameters from input cards
!        mw=vbfnloinput('#WMASS')
!        mz=vbfnloinput('#ZMASS')
!        mt=vbfnloinput('#TOPMASS')
!        mb=vbfnloinput('#BOTTOMMASS')
!        w_width=vbfnloinput('#WWIDTH')
!        z_width=vbfnloinput('#ZWIDTH')
!        nflav=vbfnloinput('#NFLAVOUR')
!        mh=vbfnloinput('#HMASS')
!        hwidth = vbfnloinput('#HWIDTH')
!        order_min=powheginput('#order_min')
!        order_max=powheginput('#order_max')
!        ebeam1=powheginput('#ebeam1')
!        ebeam2=powheginput("#ebeam2")
!        higgs_use_BW=.false.
!        if(powheginput('#higgsbreitwigner').eq.1) then
!           higgs_use_BW=.true.
!        endif
!        hmasswindow = 30.0_dp
!        if(powheginput('#higgsmasswindow').gt.0d0) then
!           hmasswindow = powheginput('#higgsmasswindow')
!        endif
!        sqrts=ebeam1+ebeam2
!        scale_choice=powheginput('#scale_choice')
!        readin=.false.
!        writeout=.false.
!        if(powheginput('#readingrid').eq.1) then
!           readin=.true.
!        endif
!        if(powheginput('#writeoutgrid').eq.1) then
!           writeout=.true.
!        endif
!        xmuf=powheginput('#facscfact')
!        xmur=powheginput('#renscfact')
!        ipdf=powheginput('#lhans1')
!        ncall1=powheginput('#ncall1')
!        ncall2=powheginput('#ncall2') 
!        itmx1=powheginput('#itmx1')
!        itmx2=powheginput('#itmx2')
!        iseed=powheginput('#iseed')
!        toypdf=.false.
!        if(powheginput('#toypdf').eq.1) then
!           toypdf=.true.
!        endif
!        toyas=powheginput('#toyalphas')
!        test_Q0=powheginput('#test_Q0')
!        muR_PDF=powheginput('#muR_PDF')
!     else
!        ! if read_from_card is false, read parameters from command line arguments
!        write(0,'(a)') '!!!! WARNING !!!!'
!        write(0,'(a)') '>> Input cards ignored, reading parameters from command line instead <<'
!        mw           = dble_val_opt("-mw",80.398_dp)
!        mz           = dble_val_opt("-mz",91.187_dp)
!        mt           = dble_val_opt("-mt",172.4_dp)
!        mb           = dble_val_opt("-mb",4.75_dp)
!        nflav        = int_val_opt ("-nf",5)
!        mh           = dble_val_opt("-mh",125.0_dp)
!        w_width      = dble_val_opt("-wwidth",2.141_dp)
!        z_width      = dble_val_opt("-zwidth",2.4952_dp)
!        hwidth       = dble_val_opt("-hwidth",0.004029643852284941_dp)
!        order_min    = int_val_opt ('-order-min',1)
!        order_max    = int_val_opt ('-order-max',3)
!        sqrts        = dble_val_opt("-sqrts",13000.0_dp)
!        scale_choice = int_val_opt ('-scale-choice',1)
!        readin       = log_val_opt ("-readingrid")
!        writeout     = log_val_opt ("-writeoutgrid")
!        higgs_use_BW = log_val_opt ("-higgsbreitwigner")
!        hmasswindow  = dble_val_opt('-higgsmasswindow',30.0_dp)
!        xmuf         = dble_val_opt("-xmuf",1.0_dp)
!        xmur         = dble_val_opt("-xmur",1.0_dp)
!        ipdf         = int_val_opt ("-pdf",261000)
!        ncall1       = int_val_opt ("-ncall1",100000)
!        ncall2       = int_val_opt ("-ncall2",100000)
!        itmx1        = int_val_opt ("-itmx1",12)
!        itmx2        = int_val_opt ("-itmx2",12)
!        iseed        = int_val_opt ("-iseed",10)
!        toypdf       = log_val_opt("-toy") ! for debugging only
!        toyas        = dble_val_opt("-toy-alphas") ! for debugging only
!        test_Q0      = dble_val_opt("-test_Q0",-1.0_dp) ! for debugging only
!        muR_PDF      = dble_val_opt("-muR_PDF") ! for debugging only
!     endif
!     
!     ! common setup
!     write(seedstr,"(I4.4)") iseed
!     idum = -iseed
!     Qmin = 1.0_dp
!     pi   = 4.0_dp*atan(1.0_dp)
!     Q0_cut_sq = 4.0_dp
!     mh_sq = mh**2
!     S = sqrts**2
!     ! compute sin(\theta_w) from W/Z mass
!     sin_thw = 1.0_dp - (mw/mz)**2
!   
!   end subroutine set_parameters

end module parameters
