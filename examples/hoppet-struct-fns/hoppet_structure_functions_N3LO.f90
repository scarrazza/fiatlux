!======================================================================
! List of issues relative to 1109.3717
! ------------------------------------
!
! - product of sums in Eq.(3.19) can't be right
!
! - getting W- from q_ns(-) + q_s would give 2u + 2dbar, which
!   interacts with W+, not W-; similar issue in F3
!
! - An additional factor of two seems to be needed in Eqs. (3.8),
!   (3.9), (3.17) and (3.18) in order get agreement with Marco
!   Zaro's code. (There he puts it in front of the vi's and ai's)
!
! - Just below Eq.(3.18) they say C_{3,ns}^V = C_{i,ns}^-: should
!   the "i" in the second term actually be a "3"?
!======================================================================
!
! Summary of our understanding of the coefficient functions
! ---------------------------------------------------------
!
! - F1 = (F2 - FL)/2x
!
! - hep-ph/0504042 (MMV): gives F2 and FL results (electromagnetic)
!
!   (1/x) F_a = C_{a,ns} \otimes q_{ns}
!                   +  <e^2> (C_{a,q} \otimes q_s  + C_{a,g} \otimes g)
!
!   with z = 2,L, where:
!
!   * <e^2> is the average electric charge
!   * q_{ns} = ???? (we're supposed to deduce it from Eq.(4.2))
!   * C_{a,q} = C_{a,ns} + C_{a,ps};
!   * C_{2,ps}^{(0)} = C_{2,ps}^{(1)} = 0 (and presumably for FL too)
!
! - http://www.sciencedirect.com/science/article/pii/055032139290087R#
!   (Zijlstra & van Neerven) has the second-order coefficient
!   functions. 
!
! - from 1109.3717 (Zaro & Co):
!
!   * q_{ns,i}^+ = (q_i + qbar_i) - q_s
!   * q_{ns,i}^- = (q_i - qbar_i) - q_{ns}^v
!   * q_{ns}^v   = \sum_{i=1}^{n_f} (q_i - qbar_i)
!   * q_s        = \sum_{i=1}^{n_f} (q_i + qbar_i)
!
!   That, together with [C_{a,q} = C_{a,ns} + C_{a,ps}] means that the combination
! 
!       q_{ns,j}^+ * C_{i,ns}^+ + q_s * C_{i,q}
!
!   reduces to 
!
!       (q_j+qbar_j) * C_{i,ns}^+ + q_s * C_{i,ps}
!
module hoppet_structure_functions
  use hoppet_v1
  use hoppet_coefficient_functions
  use parameters
  implicit none

  private
  public :: hoppet_start, hoppetStartExtendedLocal
  public :: hoppet_read_PDF
  !public :: hoppet_matrix_element
  public :: hoppet_set_structure_functions_up_to
  ! public :: hoppet_write_f1
  ! public :: hoppet_write_f2
  ! public :: hoppet_write_f3

  public :: F_LO, F_NLO, F_NNLO, F_N3LO

    real(dp), external :: alphasPDF
  type(running_coupling), save :: coupling
  !! holds information about the grid
  type(grid_def), public, save :: grid, gdarray(4)
  type(grid_def),         save :: grid_n3lo, gdarray_n3lo(4)

  public :: coupling
  
  !! holds the splitting functions
  type(dglap_holder), save :: dh

  !! 0 is main pdf table, while i=1:8 contain convolutions with the
  !! splitting function
  type(pdf_table), save :: tables(0:15)
  ! indices for the different structure functions
  integer, parameter, public :: FLWp= 1, F2Wp= 2, F3Wp= 3
  integer, parameter, public :: FLWm=-1, F2Wm=-2, F3Wm=-3
  integer, parameter, public :: FLZ = 4, F2Z = 5, F3Z = 6
  integer, parameter, public :: FLEM = -4, F2EM = -5
  
  ! constants and fixed parameters
  real(dp), parameter :: viW = 1/sqrt(two), aiW = viW
  real(dp), save      :: vi2_ai2_avg_W, two_vi_ai_avg_W
  real(dp), save      :: vi2_ai2_Z_down, vi2_ai2_Z_up
  real(dp), save      :: two_vi_ai_Z_down, two_vi_ai_Z_up
  real(dp), parameter :: e2_up = 4.0_dp/9.0_dp, e2_down = 1.0_dp/9.0_dp
  ! these log terms are only used for scale_choice = 0,1, to speed up the code
  real(dp), save      :: log_muF2_over_Q2, log_muR2_over_Q2
  
  interface operator(.dot.)
     module procedure dot
  end interface operator(.dot.)

  ! real(dp), save        :: C2LO,   CLLO,   C3LO
  ! type(split_mat), save :: C2NLO,  CLNLO,  C3NLO
  ! type(split_mat), save :: C2NNLO, CLNNLO, C3NNLO
  ! type(split_mat), save :: C2N3LO, CLN3LO, C3N3LO
  ! type(split_mat), save :: C2N3LO_fl11, CLN3LO_fl11
  ! integer :: nf_lcl

  type coeff_fn_holder
    real(dp)        :: C2LO,   CLLO,   C3LO
    type(split_mat) :: C2NLO,  CLNLO,  C3NLO
    type(split_mat) :: C2NNLO, CLNNLO, C3NNLO
    type(split_mat) :: C2N3LO, CLN3LO, C3N3LO
    type(split_mat) :: C2N3LO_fl11, CLN3LO_fl11
    integer :: nf_lcl
  end type coeff_fn_holder

  ! if this is true, then we use a zero-mass scheme that assumes
  ! a fixed number of flavours (nflav from parameters.f90)
  !
  ! Otherwise, use a zero-mass variable-flavour number scheme
  logical :: zm_ffns = .false.
  
  integer, parameter :: cfh_nflo = 3, cfh_nfhi = 6
  type(coeff_fn_holder), target,save  :: cfh_array(cfh_nflo:cfh_nfhi)
  type(coeff_fn_holder), pointer :: cfh
  

  
contains

  !subroutine hoppet_start()
  !end subroutine hoppet_start
  
  
  !----------------------------------------------------------------------
  ! This routine assumes that LHAPDF has already been initialised with a PDF
  subroutine hoppet_start()
    !real(dp), intent(in) :: 
    !----------------------------------------------------------------------
    real(dp)  :: ymax, dy, minQval, maxQval, dlnlnQ
    integer   :: nloop, order, nf_loop
    ! evaluate parameters needed for the structure functions
    ! cf. Eqs. (3.10+11) of 1109.3717
    vi2_ai2_Z_up     = one/four + (half - (four/three) * sin_thw)**2
    vi2_ai2_Z_down   = one/four + (half - (two/three)  * sin_thw)**2
    two_vi_ai_Z_up   = half - (four/three) * sin_thw
    two_vi_ai_Z_down = half - (two/three)  * sin_thw

    ! cf. Eq.3.20 + 3.17 & 3.18
    ! 1/nf \sum_j=1^nf (vj^2 + aj^2)
    vi2_ai2_avg_W = viW**2 + aiW**2
    ! 1/nf \sum_j=1^nf 2*vj*aj
    two_vi_ai_avg_W = two * viW * aiW

    ! Streamlined initialization
    ! including  parameters for x-grid
    order = -6 
    ymax  = 20.0_dp
    dy    = 0.05_dp  ! dble_val_opt("-dy",0.1_dp)
    dlnlnQ = dy/4.0_dp
    nloop = 3
    minQval = min(xmuf*Qmin, Qmin)
    maxQval = 1e6_dp
    
    call hoppetStartExtendedLocal(ymax,dy,minQval,maxQval,dlnlnQ,nloop,&
         &         order,factscheme_MSbar)

    do nf_loop = cfh_nflo, cfh_nfhi
      call SetNfCoeffFnHolder(nf_loop)
      cfh%nf_lcl = nf_int
      write(0,*) "cfh%nf_lcl = ", cfh%nf_lcl
      
      ! first initialise the LO coefficient "functions" (just a number, since delta-fn)
      cfh%C2LO = one
      cfh%CLLO = zero
      cfh%C3LO = one
      
      ! now initialise some NLO coefficient functions
      call InitC2NLO(grid, cfh%C2NLO)
      call InitCLNLO(grid, cfh%CLNLO)
      call InitC3NLO(grid, cfh%C3NLO) 
      
      ! and the NNLO ones
      call InitC2NNLO(grid, cfh%C2NNLO)
      call InitCLNNLO(grid, cfh%CLNNLO)
      call InitC3NNLO(grid, cfh%C3NNLO) 
      
      ! and the N3LO ones
      call InitC2N3LO(grid_n3lo, cfh%C2N3LO)
      call InitCLN3LO(grid_n3lo, cfh%CLN3LO)
      call InitC3N3LO(grid_n3lo, cfh%C3N3LO)
      ! including the fl11 terms for Z/photon exchanges
      call InitC2N3LO_fl11(grid_n3lo, cfh%C2N3LO_fl11)
      call InitCLN3LO_fl11(grid_n3lo, cfh%CLN3LO_fl11)
    end do
   
      
    !-- then go to the standard one we had
    call SetNfCoeffFnHolder(nflav)
  end subroutine hoppet_start

  subroutine SetNfCoeffFnHolder(nf_here)
    integer, intent(in) :: nf_here
    call SetNfDglapHolder(dh, nf_here)
    cfh => cfh_array(nf_here)
  end subroutine SetNfCoeffFnHolder
  
  
  !----------------------------------------------------------------------
  subroutine hoppetStartExtendedLocal(ymax,dy,valQmin,valQmax,dlnlnQ,nloop,order,factscheme)
    implicit none
    real(dp), intent(in) :: ymax   !! highest value of ln1/x user wants to access
    real(dp), intent(in) :: dy     !! internal grid spacing: 0.1 is a sensible value
    real(dp), intent(in) :: valQmin, valQmax !! range in Q
    real(dp), intent(in) :: dlnlnQ !! internal table spacing in lnlnQ
    integer,  intent(in) :: nloop  !! the maximum number of loops we'll want (<=3)
    integer,  intent(in) :: order  !! order of numerical interpolation (+ve v. -ve: see below)
    integer,  intent(in) :: factscheme !! 1=unpol-MSbar, 2=unpol-DIS, 3=Pol-MSbar
    !-------------------------------------

    ! initialise our grids
    if (order >= 0) call wae_error("hoppetStartExtendedLocal", "order >= 0 not allowed, order was",&
         & intval=order)
    ! the internal interpolation order (with a minus sign allows
    ! interpolation to take fake zero points beyond x=1 -- convolution
    ! times are unchanged, initialisation time is much reduced and
    ! accuracy is slightly reduced)
    !order = -5
    ! Now create a nested grid
    call InitGridDef(gdarray(4),dy/27.0_dp,0.2_dp, order=order)
    call InitGridDef(gdarray(3),dy/9.0_dp,0.5_dp,  order=order)
    call InitGridDef(gdarray(2),dy/3.0_dp,2.0_dp,  order=order)
    call InitGridDef(gdarray(1),dy,       ymax  ,  order=order)
    call InitGridDef(grid,gdarray(1:4),locked=.true.)

    ! At N3LO we need to change the grid slightly:
    !  - convolution epsilon is reduced from 10^-7 to 10^-6
    !  - interpolation order is reduced from -6 to -5 for the high x region
    ! Note: The grid points in y have to be the same as for the other grid!
    call SetDefaultConvolutionEps(0.000001_dp) ! anything less breaks integration of N3LO pieces
    call InitGridDef(gdarray_n3lo(4),dy/27.0_dp,0.2_dp, order=-abs(order)+1)!reduce order to -5
    call InitGridDef(gdarray_n3lo(3),dy/9.0_dp,0.5_dp,  order=-abs(order)+1) !reduce order to -5
    call InitGridDef(gdarray_n3lo(2),dy/3.0_dp,2.0_dp,  order=order)
    call InitGridDef(gdarray_n3lo(1),dy,       ymax  ,  order=order)
    call InitGridDef(grid_n3lo,gdarray_n3lo(1:4),locked=.true.)

    ! create the tables that will contain our copy of the user's pdf
    ! as well as the convolutions with the pdf.
    call AllocPdfTable(grid, tables(:), valQmin, valQmax, &
         & dlnlnQ = dlnlnQ, freeze_at_Qmin=.true.)
    
    ! tables(0) : pdf
    tables(0)%tab = zero
    ! tables(1) : LO structure function : C_LO x f
    tables(1)%tab = zero
    
    ! tables(2) : NLO structure function : C_NLO x f
    tables(2)%tab = zero
    ! tables(3) : NLO contribution : C_LO x P_LO x f
    tables(3)%tab = zero
    
    ! tables(4) : NNLO structure function : C_NNLO x f
    tables(4)%tab = zero
    ! tables(5) : NNLO contribution : C_LO x P_LO^2 x f
    tables(5)%tab = zero
    ! tables(6) : NNLO contribution : C_LO x P_NLO x f
    tables(6)%tab = zero
    ! tables(7) : NNLO contribution : C_NLO x P_LO x f
    tables(7)%tab = zero
    
    ! tables(8)  : N3LO contribution : C_N3LO x f
    tables(8)%tab = zero
    ! tables(9)  : N3LO contribution : C_LO x P_LO^3 x f
    tables(9)%tab = zero
    ! tables(10) : N3LO contribution : C_LO x P_LO x P_NLO x f
    tables(10)%tab = zero
    ! tables(11) : N3LO contribution : C_LO x P_NLO x P_LO x f
    tables(11)%tab = zero
    ! tables(12) : N3LO contribution : C_NLO x P_LO^2 x f
    tables(12)%tab = zero
    ! tables(13) : N3LO contribution : C_NLO x P_NLO x f
    tables(13)%tab = zero
    ! tables(14) : N3LO contribution : C_NNLO x P_LO x f
    tables(14)%tab = zero
    ! tables(15) : N3LO contribution : C_LO x P_NNLO x f
    tables(15)%tab = zero
    
    ! initialise splitting-function holder
    call InitDglapHolder(grid,dh,factscheme=factscheme,&
         &                      nloop=nloop,nflo=3,nfhi=6)
    ! choose a sensible default number of flavours.
    call SetNfDglapHolder(dh,nflcl=nflav)
  end subroutine hoppetStartExtendedLocal


  !----------------------------------------------------------------------
  ! fill the streamlined interface PDF table (possibly using hoppet's
  ! evolution)
  subroutine hoppet_read_PDF()
    use helpers
    !real(dp) :: muR_Q
    !real(dp) :: Q0pdf, asMZ
    !! define the interfaces for LHA pdf (by default not used)
    !! (NB: unfortunately this conflicts with an internal hoppet name,
    !! so make sure that you "redefine" the internal hoppet name,
    !! as illustrated in the commented "use" line above:
    !! use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, ...)
    interface
       subroutine EvolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine EvolvePDF
    end interface
    integer :: i
    !logical, save :: pre_ev_done = .false.
    !----------------
    real(dp) :: res_lhapdf(-6:6), x, Q
    real(dp) :: res_hoppet(-6:6)
    real(dp) :: pdf_at_Q0(0:grid%ny,ncompmin:ncompmax)
    integer  :: ix
    real(dp) :: xvals(0:grid%ny), quark_masses(4:6)
    
    ! ! Set parameters of running coupling
    ! asMZ = alphasPDF(MZ)
    ! Q0pdf = 10.0_dp
    ! muR_Q = 1.0_dp
    ! if (hoppet_evolution) then
    !    if (hoppet_pre_evolution) then
    !       if (.not. pre_ev_done) then
    !          call hoppetPreEvolve(asMZ, MZ, nloop, muR_Q, Q0pdf)
    !          pre_ev_done = .true.
    !       end if
    !       call hoppetCachedEvolve(EvolvePDF)
    !    else
    !       call hoppetEvolve(asMZ, MZ, nloop, muR_Q, EvolvePDF, Q0pdf)
    !    end if
    !    write(6,'(a)') "Evolution done!"
    ! else
    !    call hoppetAssign(EvolvePDF)
    !    write(6,'(a)') "PDF assigned to hoppet table!"
    ! end if
    ! !write(0,*) "Size of table = ", size(tables(0)%tab)

    if (test_Q0 > zero) then

       write(0,*) "WARNING: Using internal HOPPET DGLAP evolution"
       xvals = xValues(grid)

       do ix = 0, grid%ny
          call EvolvePDF(xvals(ix),test_Q0,pdf_at_Q0(ix,:))
       enddo
       call InitRunningCoupling(coupling, alphasPDF(MZ) , MZ , 4,&
            & -1000000045, (/  1.275_dp, 4.18_dp, 173.21_dp/),&
            & .true.)
       call EvolvePdfTable(tables(0), test_Q0, pdf_at_Q0, dh, coupling, &
            &  muR_Q=muR_PDF, nloop=3)
    else
    ! InitRunningCoupling has to be called for the HOPPET coupling to be initialised 
    ! Default is to ask for 4 loop running and threshold corrections at quark masses.  
       quark_masses(4:6) = (/lhapdf_qmass(4), lhapdf_qmass(5), lhapdf_qmass(6)/)
       ! manually override top mass, since lhapdf usually doesn't have 6 flavours, only 5
       quark_masses(6) = 1e100_dp
       write(6,*) 'quark masses are = ', quark_masses
       call InitRunningCoupling(coupling, alfas=alphasPDF(MZ), Q=MZ ,nloop=4,&
            & quark_masses = quark_masses(4:6),&
            & masses_are_MSbar = .false., muMatch_mQuark = muR_PDF)

       ! now fix up nf info
       do i = lbound(tables,1), ubound(tables,1)
          call AddNfInfoToPdfTable(tables(i), coupling)
          tables(i)%tab = zero
       end do
       
       call hoppetAssign(EvolvePDF)
    end if


    
    ! quickly test that we have read in the PDFs correctly
    write(0,*) "Quick test that PDFs have been read in correctly"
    x = 0.08_dp
    Q = 17.0_dp
    call EvolvePDF(x, Q, res_lhapdf)
    write(0,*) 'Test that LHAPDF and HOPPET PDFs agree'
    write(0,*) '  lhapdf: ', x, Q, res_lhapdf
    call EvalPdfTable_xQ(tables(0), x, Q, res_hoppet)
    write(0,*) '  hoppet: ', x, Q, res_hoppet
  end subroutine hoppet_read_PDF

  !----------------------------------------------------------------------
  !! Given a pdf_subroutine with the interface shown below, initialise
  !! our internal pdf table.
  subroutine hoppetAssign(pdf_subroutine)
    implicit none
    interface ! indicate what "interface" pdf_subroutine is expected to have
       subroutine pdf_subroutine(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine pdf_subroutine
    end interface
    !-----------------------------------

    ! set up table(0) by copying the values returned by pdf_subroutine onto 
    ! the x-Q grid in table(0)
    call FillPdfTable_LHAPDF(tables(0), pdf_subroutine)
  end subroutine hoppetAssign

  !----------------------------------------------------------------------
  ! Set up the structure function tables up to a given order
  subroutine hoppet_set_structure_functions_up_to(order)
    integer, intent(in) :: order

    ! To turn off b quarks completely (only for testing and comparison)
    ! uncomment following two lines:
    ! tables(0)%tab(:,-5,:) = zero
    ! tables(0)%tab(:,+5,:) = zero
    
    !if (scale_choice.le.1) then
       ! if scale_choice = 0,1 use fast implementation
       ! tables is saved as an array in Q, and only tables(0),
       ! tables(1), tables(2), tables(4), tables(8) are non zero
       if (order.ge.1) call hoppet_set_LO_structure_functions()
       if (order.ge.2) call hoppet_set_NLO_structure_functions()
       if (order.ge.3) call hoppet_set_NNLO_structure_functions()
       if (order.ge.4) call hoppet_set_N3LO_structure_functions()
    !else
    !   ! if scale_choice >= 2 use slower implementation with full
    !   ! scale choices, such as sqrt(Q1*Q2)
    !   ! tables is saved as an array in muF now, and all components
    !   ! of tables are non zero.
    !   if (order.ge.1) call hoppet_set_LO_structure_functions_anyscale()
    !   if (order.ge.2) call hoppet_set_NLO_structure_functions_anyscale()
    !   if (order.ge.3) call hoppet_set_NNLO_structure_functions_anyscale()
    !   if (order.ge.4) call hoppet_set_N3LO_structure_functions_anyscale()
    !endif

    ! To rescale PDF by the N3LO F2 structure function (evaluated at 10 GeV)
    ! as a check of the size of missing N3LO PDFs, uncomment following line:
    ! call rescale_pdf_nnlo(10.0_dp)
    ! call rescale_pdf_n3lo(10.0_dp)
  end subroutine hoppet_set_structure_functions_up_to

!   !----------------------------------------------------------------------
!   ! Rescale the PDF by F2 NNLO / F2 N3LO
!   subroutine rescale_pdf_n3lo (Qval)
!     real(dp), intent(in) :: Qval
!     real(dp) :: muF, muR, factor
!     real(dp) :: str_fct(-6:6), f2nnlo, f2n3lo, y(0:grid%ny)
!     integer :: iy, ii
!     muF = muF1(Qval,Qval, 0.0_dp)
!     muR = muR1(Qval,Qval, 0.0_dp)
!     y = yValues(grid)
!     do iy = 0, grid%ny
!        str_fct(:) = F_LO(y(iy), Qval, muR, muF) + F_NLO(y(iy), Qval, muR, muF) + F_NNLO(y(iy), Qval, muR, muF)
!        f2nnlo = str_fct(F2Z)
!        str_fct(:) = str_fct(:) + F_N3LO(y(iy), Qval, muR, muF)
!        f2n3lo = str_fct(F2Z)
!        factor = 1.0_dp
!        if (f2n3lo.gt.0.0_dp) factor = f2nnlo/f2n3lo
!        do ii = ncompmin, ncompmax
!           tables(0)%tab(iy,ii,:) = tables(0)%tab(iy,ii,:) * factor
!        enddo
!        
!     enddo
! 
!     ! Re-set the structure functions using updated PDF
!     if (scale_choice.le.1) then
!        call hoppet_set_LO_structure_functions()
!        call hoppet_set_NLO_structure_functions()
!        call hoppet_set_NNLO_structure_functions()
!        call hoppet_set_N3LO_structure_functions()
!     else
!        call hoppet_set_LO_structure_functions_anyscale()
!        call hoppet_set_NLO_structure_functions_anyscale()
!        call hoppet_set_NNLO_structure_functions_anyscale()
!        call hoppet_set_N3LO_structure_functions_anyscale()
!     endif
! 
!   end subroutine rescale_pdf_n3lo
  
!  !---------------------------------------------------------------------- 
!  ! Rescale the PDF by F2 NLO / F2 NNLO
!  subroutine rescale_pdf_nnlo (Qval)
!    real(dp), intent(in) :: Qval
!    real(dp) :: muF, muR, factor
!    real(dp) :: str_fct(-6:6), f2nlo, f2nnlo, y(0:grid%ny) 
!    integer :: iy, ii
!    muF = muF1(Qval,Qval, 0.0_dp) 
!    muR = muR1(Qval,Qval, 0.0_dp)
!    y = yValues(grid)
!    do iy = 0, grid%ny              
!       str_fct(:) = F_LO(y(iy), Qval, muR, muF) + F_NLO(y(iy), Qval, muR, muF)
!       f2nlo = str_fct(F2Z) 
!       str_fct(:) = str_fct(:) + F_NNLO(y(iy), Qval, muR, muF)
!       f2nnlo = str_fct(F2Z)
!       factor = 1.0_dp
!       if (f2nnlo.gt.0.0_dp) factor = f2nlo/f2nnlo
!       do ii = ncompmin, ncompmax
!          tables(0)%tab(iy,ii,:) = tables(0)%tab(iy,ii,:) * factor
!       enddo
!    enddo
!    
!    ! Re-set the structure functions using updated PDF
!    if (scale_choice.le.1) then
!       call hoppet_set_LO_structure_functions()
!       call hoppet_set_NLO_structure_functions()
!       call hoppet_set_NNLO_structure_functions()
!       call hoppet_set_N3LO_structure_functions()
!    else
!       call hoppet_set_LO_structure_functions_anyscale()
!       call hoppet_set_NLO_structure_functions_anyscale()
!       call hoppet_set_NNLO_structure_functions_anyscale()
!       call hoppet_set_N3LO_structure_functions_anyscale()
!    endif
!  end subroutine rescale_pdf_nnlo
  
!   !----------------------------------------------------------------------
!   ! write the F1 structure function to idev
!   subroutine hoppet_write_f1 (idev, Q1test, Q2test, ymax, ny)
!     real(dp), intent(in) :: Q1test, Q2test, ymax
!     integer, intent(in)  :: idev, ny
!     real(dp) :: ytest, xval, muR, muF, F1Z_LO, F1Z_NLO, F1Z_NNLO, F1Z_N3LO, res(-6:6)
!     integer  :: iy
!     !F1 Wp Wm Z
!     write(idev,'(a,f10.4,a,f10.4)') '# Q1 = ', Q1test,' , Q2 = ', Q2test
!     write(idev,'(a,a)') '# x  F1Wp(LO) F1Wm(LO) F1Wp(NLO) F1Wm(NLO) F1Wp(NNLO) F1Wm(NNLO)', &
!          & ' F1Wp(N3LO) F1Wm(N3LO) F1Z(LO) F1Z(NLO) F1Z(NNLO) F1Z(N3LO)'
!     muF = muF1(Q1test,Q2test, 0.0_dp)
!     muR = muR1(Q1test,Q2test, 0.0_dp)
!     do iy = ny, 1, -1
!        ytest = iy * ymax / ny
!        xval = exp(-ytest)
!        res = F_LO(ytest, Q1test, muR, muF)
!        write(idev,'(3es22.12)',advance='no') xval, res(F1Wp),res(F1Wm)
!        F1Z_LO = res(F1Z)
!        res = F_NLO(ytest, Q1test, muR, muF)
!        write(idev,'(2es22.12)',advance='no') res(F1Wp), res(F1Wm)
!        F1Z_NLO = res(F1Z)
!        res = F_NNLO(ytest, Q1test, muR, muF)
!        write(idev,'(2es22.12)',advance='no') res(F1Wp), res(F1Wm)
!        F1Z_NNLO = res(F1Z)
!        res = F_N3LO(ytest, Q1test, muR, muF)
!        write(idev,'(2es22.12)',advance='no') res(F1Wp), res(F1Wm)
!        F1Z_N3LO = res(F1Z)
!        write(idev,'(4es22.12)',advance='no') F1Z_LO, F1Z_NLO, F1Z_NNLO, F1Z_N3LO
!        write(idev,*)
!     end do
!   end subroutine hoppet_write_f1
!   
!   !----------------------------------------------------------------------
!   ! write the F2 structure function to idev
!   subroutine hoppet_write_f2 (idev, Q1test, Q2test, ymax, ny)
!     real(dp), intent(in) :: Q1test, Q2test, ymax
!     integer, intent(in)  :: idev, ny
!     real(dp) :: ytest, xval, muR, muF, F2Z_LO, F2Z_NLO, F2Z_NNLO, F2Z_N3LO, res(-6:6)
!     integer  :: iy
!     !F2 Wp Wm Z
!     write(idev,'(a,f10.4,a,f10.4)') '# Q1 = ', Q1test,' , Q2 = ', Q2test
!     write(idev,'(a,a)') '# x  F2Wp(LO) F2Wm(LO) F2Wp(NLO) F2Wm(NLO) F2Wp(NNLO) F2Wm(NNLO)', &
!          & ' F2Wp(N3LO) F2Wm(N3LO) F2Z(LO) F2Z(NLO) F2Z(NNLO) F2Z(N3LO)'
!     muF = muF1(Q1test,Q2test, 0.0_dp)
!     muR = muR1(Q1test,Q2test, 0.0_dp)
!     do iy = ny, 1, -1
!        ytest = iy * ymax / ny
!        xval = exp(-ytest)
!        res = F_LO(ytest, Q1test, muR, muF)
!        write(idev,'(3es22.12)',advance='no') xval, res(F2Wp),res(F2Wm)
!        F2Z_LO = res(F2Z)
!        res = F_NLO(ytest, Q1test, muR, muF)
!        write(idev,'(2es22.12)',advance='no') res(F2Wp), res(F2Wm)
!        F2Z_NLO = res(F2Z)
!        res = F_NNLO(ytest, Q1test, muR, muF)
!        write(idev,'(2es22.12)',advance='no') res(F2Wp), res(F2Wm)
!        F2Z_NNLO = res(F2Z)
!        res = F_N3LO(ytest, Q1test, muR, muF)
!        write(idev,'(2es22.12)',advance='no') res(F2Wp), res(F2Wm)
!        F2Z_N3LO = res(F2Z)
!        write(idev,'(4es22.12)',advance='no') F2Z_LO, F2Z_NLO, F2Z_NNLO, F2Z_N3LO
!        write(idev,*)
!     end do
!   end subroutine hoppet_write_f2
! 
!   !----------------------------------------------------------------------
!   ! write the F3 structure function to idev
!   subroutine hoppet_write_f3 (idev, Q1test, Q2test, ymax, ny)
!     real(dp), intent(in) :: Q1test, Q2test, ymax
!     integer, intent(in)  :: idev, ny
!     real(dp) :: ytest, xval, muR, muF, F3Z_LO, F3Z_NLO, F3Z_NNLO, F3Z_N3LO, res(-6:6)
!     integer  :: iy
!     !F3 Wp Wm Z
!     write(idev,'(a,f10.4,a,f10.4)') '# Q1 = ', Q1test,' , Q2 = ', Q2test
!     write(idev,'(a,a)') '# x  F3Wp(LO) F3Wm(LO) F3Wp(NLO) F3Wm(NLO) F3Wp(NNLO) F3Wm(NNLO)', &
!          & ' F3Wp(N3LO) F3Wm(N3LO) F3Z(LO) F3Z(NLO) F3Z(NNLO) F3Z(N3LO)'
!     muF = muF1(Q1test,Q2test, 0.0_dp)
!     muR = muR1(Q1test,Q2test, 0.0_dp)
!     do iy = ny, 1, -1
!        ytest = iy * ymax / ny
!        xval = exp(-ytest)
!        res = F_LO(ytest, Q1test, muR, muF)
!        write(idev,'(3es22.12)',advance='no') xval, res(F3Wp),res(F3Wm)
!        F3Z_LO = res(F3Z)
!        res = F_NLO(ytest, Q1test, muR, muF)
!        write(idev,'(2es22.12)',advance='no') res(F3Wp), res(F3Wm)
!        F3Z_NLO = res(F3Z)
!        res = F_NNLO(ytest, Q1test, muR, muF)
!        write(idev,'(2es22.12)',advance='no') res(F3Wp), res(F3Wm)
!        F3Z_NNLO = res(F3Z)
!        res = F_N3LO(ytest, Q1test, muR, muF)
!        write(idev,'(2es22.12)',advance='no') res(F3Wp), res(F3Wm)
!        F3Z_N3LO = res(F3Z)
!        write(idev,'(4es22.12)',advance='no') F3Z_LO, F3Z_NLO, F3Z_NNLO, F3Z_N3LO
!        write(idev,*)
!     end do
!   end subroutine hoppet_write_f3


  !----------------------------------------------------------------------
  ! make sure we are set up for the nf value corresponding to entry iQ
  ! of the table
  subroutine SwitchNfFromTable(iQ)
    integer,         intent(in) :: iQ
    if (.not. zm_ffns .and. tables(0)%nf_info_associated) &
            &call SetNfCoeffFnHolder(tables(0)%nf_int(iQ))
  end subroutine SwitchNfFromTable
  
  
  !----------------------------------------------------------------------
  ! set up LO structure functions for scale_choice = 0, 1
  subroutine hoppet_set_LO_structure_functions()
    integer :: iQ
    real(dp) :: f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q

    ! all coefficient functions at LO are delta functions (F2, FL and F3),
    ! so simply pass table(0) for each of the pieces
    do iQ = 0, tables(0)%nQ
       call SwitchNfFromTable(iQ)
       Q = tables(0)%Q_vals(iQ)

       ! explicitly evaluate the PDF at scale muF(Q)
       call EvalPdfTable_Q(tables(0),muF(Q),f)
       tables(1)%tab(:,:,iQ) = structure_function_general(cfh%C2LO*f, cfh%CLLO*f, cfh%C3LO*f)
    end do
    
  end subroutine hoppet_set_LO_structure_functions
  
  !----------------------------------------------------------------------
  ! Set up the LO structure functions for any scale choice
  subroutine hoppet_set_LO_structure_functions_anyscale()
    integer :: iQ

    ! all coefficient functions at LO are delta functions (F2, FL and F3),
    ! so simply pass table(0) for each of the pieces
    do iQ = 0, tables(0)%nQ
       call SwitchNfFromTable(iQ)
       tables(1)%tab(:,:,iQ) = structure_function_general(&
            & tables(0)%tab(:,:,iQ) * cfh%C2LO, &
            & tables(0)%tab(:,:,iQ) * cfh%CLLO, &
            & tables(0)%tab(:,:,iQ) * cfh%C3LO)       
    end do
    
  end subroutine hoppet_set_LO_structure_functions_anyscale
  
  !----------------------------------------------------------------------
  ! set up NLO structure functions for scale_choice = 0, 1
  subroutine hoppet_set_NLO_structure_functions()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q
    
    do iQ = 0, tables(0)%nQ

       call SwitchNfFromTable(iQ)
       
       Q = tables(0)%Q_vals(iQ)
       call EvalPdfTable_Q(tables(0),muF(Q),f)
       call set_scale_logs(Q)
       ! do the convolution with the coefficient functions and also the
       ! corresponding splitting-function contributions when scales
       ! are not equal to Q
       PLO_f = dh%P_LO * f
       f2 = CxNLO_with_logs(cfh%C2LO, cfh%C2NLO, f, PLO_f)
       fL = CxNLO_with_logs(cfh%CLLO, cfh%CLNLO, f, PLO_f)
       f3 = CxNLO_with_logs(cfh%C3LO, cfh%C3NLO, f, PLO_f)
       
       tables(2)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

    end do

  end subroutine hoppet_set_NLO_structure_functions

  !----------------------------------------------------------------------
  ! returns the convolution of coefficient and splitting functions
  ! for NLO, with scale dependence; leading factor of as2pi left out.
  !
  ! This routine assumes that set_scale_logs(Q) has been called
  ! beforehand.
  function CxNLO_with_logs(CxLO, CxNLO, f, PLO_f) result(res)
    real(dp),        intent(in) :: CxLO
    type(split_mat), intent(in) :: CxNLO
    real(dp),        intent(in) :: f    (0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp)                    :: res  (0:grid%ny,ncompmin:ncompmax)
    !----------------------------------------------------------------

    res = CxNLO * f - (CxLO * log_muF2_over_Q2) * PLO_f
  end function CxNLO_with_logs
  
  !----------------------------------------------------------------------
  ! Set up the NLO structure functions for any scale choice
  subroutine hoppet_set_NLO_structure_functions_anyscale()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q
    
    do iQ = 0, tables(0)%nQ
       call SwitchNfFromTable(iQ)

       ! Internal Q value effectively corresponds to muF(Q1,Q2)
       Q = tables(0)%Q_vals(iQ)
       f = tables(0)%tab(:,:,iQ)

       ! Save the NLO pieces in tables(2) and tables(3)
      
       ! Get the NLO coefficient function, (C_NLO x f) 
       f2 = (cfh%C2NLO * f)
       fL = (cfh%CLNLO * f)
       f3 = (cfh%C3NLO * f)
       tables(2)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now compute the (C_LO x P_LO x f) term
       PLO_f = dh%P_LO * f
       f2 = (cfh%C2LO * PLO_f)
       fL = (cfh%CLLO * PLO_f)
       f3 = (cfh%C3LO * PLO_f)
       tables(3)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)
    end do

  end subroutine hoppet_set_NLO_structure_functions_anyscale
  
  !----------------------------------------------------------------------
  ! set up the NNLO structure functions for scale_choice = 0, 1
  subroutine hoppet_set_NNLO_structure_functions()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q
    
    do iQ = 0, tables(0)%nQ
       call SwitchNfFromTable(iQ)

       Q = tables(0)%Q_vals(iQ)
       call EvalPdfTable_Q(tables(0),muF(Q),f)
       call set_scale_logs(Q)
       ! do the convolution with the coefficient functions and also the
       ! corresponding splitting-function contributions when scales
       ! are not equal to Q
       PLO_f   = dh%P_LO  * f
       PNLO_f  = dh%P_NLO * f
       PLO2_f  = dh%P_LO  * PLO_f
       f2 = CxNNLO_with_logs(cfh%C2LO, cfh%C2NLO, cfh%C2NNLO, f, PLO_f, PNLO_f, PLO2_f)
       fL = CxNNLO_with_logs(cfh%CLLO, cfh%CLNLO, cfh%CLNNLO, f, PLO_f, PNLO_f, PLO2_f)
       f3 = CxNNLO_with_logs(cfh%C3LO, cfh%C3NLO, cfh%C3NNLO, f, PLO_f, PNLO_f, PLO2_f)

       tables(4)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)
     
    end do

  end subroutine hoppet_set_NNLO_structure_functions
  
  !----------------------------------------------------------------------
  ! returns the convolution of coefficient and splitting functions
  ! for NNLO, with scale dependence; leading factor of as2pi left out.
  !
  ! This routine assumes that set_scale_logs(Q) has been called
  ! beforehand.
  function CxNNLO_with_logs(CxLO, CxNLO, CxNNLO, f, PLO_f, PNLO_f, PLO2_f) result(res)
    real(dp),        intent(in) :: CxLO
    type(split_mat), intent(in) :: CxNLO
    type(split_mat), intent(in) :: CxNNLO
    real(dp),        intent(in) :: f    (0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp)                    :: res  (0:grid%ny,ncompmin:ncompmax)
    !----------------------------------------------------------------
    real(dp) :: f_NLO(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f_NNLO(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: LR, LF

    LR = log_muR2_over_Q2
    LF = log_muF2_over_Q2
    
    f_NLO  = - LF * PLO_f

    f_NNLO = + half * LF**2 * PLO2_f &
         &   - LF * PNLO_f &
         &   - (twopi*beta0*(LR*LF - half * LF**2)) * PLO_f
    
    res = CxNNLO * f  +  CxNLO * f_NLO  +  CxLO * f_NNLO  +  (twopi*beta0*LR) * (CxNLO * f)
  end function CxNNLO_with_logs


  !----------------------------------------------------------------------
  ! set up the NNLO structure functions
  subroutine hoppet_set_NNLO_structure_functions_anyscale()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q

    do iQ = 0, tables(0)%nQ
       call SwitchNfFromTable(iQ)

       Q = tables(0)%Q_vals(iQ)
       f = tables(0)%tab(:,:,iQ)
       
       ! save the NNLO pieces in tables(4:7)
       
       PLO2_f  = dh%P_LO  * (dh%P_LO * f)
       PLO_f   = dh%P_LO  * f
       PNLO_f  = dh%P_NLO * f

       ! first calculate the pure NNLO term, (C_NNLO x f) 
       f2 = (cfh%C2NNLO * f)
       fL = (cfh%CLNNLO * f)
       f3 = (cfh%C3NNLO * f)
       tables(4)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_LO^2 x f) term
       f2 =  (cfh%C2LO * PLO2_f)
       fL =  (cfh%CLLO * PLO2_f)
       f3 =  (cfh%C3LO * PLO2_f)
       tables(5)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_NLO) term
       f2 = (cfh%C2LO * PNLO_f)
       fL = (cfh%CLLO * PNLO_f)
       f3 = (cfh%C3LO * PNLO_f)
       tables(6)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NLO x P_LO) term
       f2 = (cfh%C2NLO * PLO_f)
       fL = (cfh%CLNLO * PLO_f)
       f3 = (cfh%C3NLO * PLO_f)
       tables(7)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

    end do

  end subroutine hoppet_set_NNLO_structure_functions_anyscale

  !----------------------------------------------------------------------
  ! set up the N3LO structure functions for scale_choice = 0, 1
  ! Warning : for now factorisation and renormalisation scale are set to Q
  subroutine hoppet_set_N3LO_structure_functions()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLONLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLOLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO3_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q
    
    do iQ = 0, tables(0)%nQ
       call SwitchNfFromTable(iQ)

       Q = tables(0)%Q_vals(iQ)
       call EvalPdfTable_Q(tables(0),muF(Q),f)
       call set_scale_logs(Q)
       ! do the convolution with the coefficient functions and also the
       ! corresponding splitting-function contributions when scales
       ! are not equal to Q
       PLO_f    = dh%P_LO   * f
       PNLO_f   = dh%P_NLO  * f
       PLO2_f   = dh%P_LO   * PLO_f
       PNNLO_f  = dh%P_NNLO * f
       PLONLO_f = dh%P_LO   * PNLO_f
       PNLOLO_f = dh%P_NLO  * PLO_f
       PLO3_f   = dh%P_LO   * PLO2_f
    
       f2 = CxN3LO_with_logs(cfh%C2LO, cfh%C2NLO, cfh%C2NNLO, cfh%C2N3LO, f, PLO_f, PNLO_f, PLO2_f, &
            &                PNNLO_f, PLONLO_f, PNLOLO_f, PLO3_f)
       fL = CxN3LO_with_logs(cfh%CLLO, cfh%CLNLO, cfh%CLNNLO, cfh%CLN3LO, f, PLO_f, PNLO_f, PLO2_f, &
            &                PNNLO_f, PLONLO_f, PNLOLO_f, PLO3_f)
       f3 = CxN3LO_with_logs(cfh%C3LO, cfh%C3NLO, cfh%C3NNLO, cfh%C3N3LO, f, PLO_f, PNLO_f, PLO2_f, &
            &                PNNLO_f, PLONLO_f, PNLOLO_f, PLO3_f)
       ! now the fl_11 piece
       f2_fl11 = cfh%C2N3LO_fl11 * f
       fL_fl11 = cfh%CLN3LO_fl11 * f

       !tables(8)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)
       tables(8)%tab(:,:,iQ) = structure_function_general_full(f2, fL, f3, f2_fl11, fL_fl11)
     
    end do

  end subroutine hoppet_set_N3LO_structure_functions
  

  !----------------------------------------------------------------------
  ! returns the convolution of coefficient and splitting functions
  ! for N3LO, with scale dependence; leading factor of as2pi left out.
  !
  ! This routine assumes that set_scale_logs(Q) has been called
  ! beforehand.
  function CxN3LO_with_logs(CxLO, CxNLO, CxNNLO, CxN3LO, f, PLO_f, PNLO_f, PLO2_f, &
       &                    PNNLO_f, PLONLO_f, PNLOLO_f, PLO3_f) result(res)
    real(dp),        intent(in) :: CxLO
    type(split_mat), intent(in) :: CxNLO
    type(split_mat), intent(in) :: CxNNLO
    type(split_mat), intent(in) :: CxN3LO
    real(dp),        intent(in) :: f    (0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PNNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLONLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PNLOLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLO3_f(0:grid%ny,ncompmin:ncompmax)
    real(dp)                    :: res  (0:grid%ny,ncompmin:ncompmax)
    !----------------------------------------------------------------
    real(dp) :: f_NLO(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f_NNLO(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f_N3LO(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: LR, LF

    LR = log_muR2_over_Q2
    LF = log_muF2_over_Q2
    
    f_NLO  = - LF * PLO_f

    f_NNLO = + half * LF**2 * PLO2_f &
         &   - LF * PNLO_f &
         &   - (twopi*beta0*(LR*LF - half * LF**2)) * PLO_f

    f_N3LO = -(1.0_dp/6.0_dp) * LF * ( &
         & - three * twopi**2 * beta1 * (LF - two * LR) * PLO_f &
         & + two * (twopi*beta0)**2 * (LF**2 - three * LF * LR + three * LR**2) * PLO_f &
         & + LF**2 * PLO3_f &
         & + three * twopi * beta0 * (LF - two * LR) * (LF * PLO2_f - two * PNLO_f) &
         & - 3.0_dp * LF * (PLONLO_f + PNLOLO_f) + 6.0_dp * PNNLO_f)
    
    res = CxN3LO * f  &
         & + (1.0_dp/6.0_dp) * LR * (6.0_dp * twopi**2 * beta1 * (CxNLO * f) &
         & + twopi*beta0 * (12.0_dp * (CxNNLO * f) + 6.0_dp * twopi * beta0 * LR * (CxNLO * f))) &
         & + CxNNLO * f_NLO + twopi * beta0 * LR * (CxNLO * f_NLO) &
         & + CxNLO * f_NNLO + CxLO * f_N3LO
  end function CxN3LO_with_logs

  !----------------------------------------------------------------------
  ! set up the N3LO structure functions
  subroutine hoppet_set_N3LO_structure_functions_anyscale()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLONLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLOLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO3_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q

    do iQ = 0, tables(0)%nQ
       call SwitchNfFromTable(iQ)

       Q = tables(0)%Q_vals(iQ)
       f = tables(0)%tab(:,:,iQ)
       
       ! save the N3LO pieces in tables(8:15)
       
       PLO2_f   = dh%P_LO  * (dh%P_LO * f)
       PLO_f    = dh%P_LO  * f
       PNLO_f   = dh%P_NLO * f
       PNNLO_f  = dh%P_NNLO * f
       PLONLO_f = dh%P_LO   * PNLO_f
       PNLOLO_f = dh%P_NLO  * PLO_f
       PLO3_f   = dh%P_LO   * (dh%P_LO * PLO_f)

       ! first calculate the pure N3LO term, (C_N3LO x f) 
       f2 = (cfh%C2N3LO * f)
       fL = (cfh%CLN3LO * f)
       f3 = (cfh%C3N3LO * f)
       f2_fl11 = (cfh%C2N3LO_fl11 * f)
       fL_fl11 = (cfh%CLN3LO_fl11 * f)
       ! tables(8)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)
       tables(8)%tab(:,:,iQ) = structure_function_general_full(f2, fL, f3, f2_fl11, fL_fl11)

       ! Now calculate the (C_LO x P_LO^3 x f) term
       f2 =  (cfh%C2LO * PLO3_f)
       fL =  (cfh%CLLO * PLO3_f)
       f3 =  (cfh%C3LO * PLO3_f)
       tables(9)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_LO x P_NLO x f) term
       f2 =  (cfh%C2LO * PLONLO_f)
       fL =  (cfh%CLLO * PLONLO_f)
       f3 =  (cfh%C3LO * PLONLO_f)
       tables(10)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_NLO x P_LO x f) term
       f2 =  (cfh%C2LO * PNLOLO_f)
       fL =  (cfh%CLLO * PNLOLO_f)
       f3 =  (cfh%C3LO * PNLOLO_f)
       tables(11)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NLO x P_LO^2 x f) term
       f2 =  (cfh%C2NLO * PLO2_f)
       fL =  (cfh%CLNLO * PLO2_f)
       f3 =  (cfh%C3NLO * PLO2_f)
       tables(12)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NLO x P_NLO) term
       f2 = (cfh%C2NLO * PNLO_f)
       fL = (cfh%CLNLO * PNLO_f)
       f3 = (cfh%C3NLO * PNLO_f)
       tables(13)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_NNLO x P_LO) term
       f2 = (cfh%C2NNLO * PLO_f)
       fL = (cfh%CLNNLO * PLO_f)
       f3 = (cfh%C3NNLO * PLO_f)
       tables(14)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

       ! Now calculate the (C_LO x P_NNLO) term
       f2 = (cfh%C2LO * PNNLO_f)
       fL = (cfh%CLLO * PNNLO_f)
       f3 = (cfh%C3LO * PNNLO_f)
       tables(15)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

    end do

  end subroutine hoppet_set_N3LO_structure_functions_anyscale


 !----------------------------------------------------------------------
 ! Structure function valid up to NNLO
  function structure_function_general(C2_f, CL_f, C3_f) result(res)
    real(dp), intent(in) :: C2_f(0:,ncompmin:), CL_f(0:,ncompmin:), C3_f(0:,ncompmin:)
    real(dp) :: C2_f_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: CL_f_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: res(0:ubound(C2_f,dim=1), ncompmin:ncompmax)
    C2_f_fl11 = zero
    CL_f_fl11 = zero
    res = structure_function_general_full(C2_f, CL_f, C3_f, C2_f_fl11, CL_f_fl11)
  end function structure_function_general


 !----------------------------------------------------------------------
 ! Structure function including N3LO fl_11 piece for neutral current
  function structure_function_general_full(C2_f, CL_f, C3_f, C2_f_fl11, CL_f_fl11) result(res)
    real(dp), intent(in) :: C2_f(0:,ncompmin:), CL_f(0:,ncompmin:), C3_f(0:,ncompmin:)
    real(dp), intent(in) :: C2_f_fl11(0:,ncompmin:), CL_f_fl11(0:,ncompmin:)
    real(dp) :: C2_f_NC(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: CL_f_NC(0:grid%ny,ncompmin:ncompmax)
    real(dp)             :: res(0:ubound(C2_f,dim=1), ncompmin:ncompmax)
    !----------------------
    ! not just up and down, but the sum 
    real(dp) :: u(0:grid%ny), d(0:grid%ny), ubar(0:grid%ny), dbar(0:grid%ny) 
    real(dp) :: two_xvals(0:grid%ny)
    integer  :: nf_save
    
    two_xvals = two*xValues(grid) ! move this outside at some point
    C2_f_NC = C2_f + C2_f_fl11
    CL_f_NC = CL_f + CL_f_fl11

    res = zero
    !--- deal with Z case -----------------------------------------
    res(:,F2Z) = (dlike(C2_f_NC) + dbarlike(C2_f_NC))*vi2_ai2_Z_down + &
         &       (ulike(C2_f_NC) + ubarlike(C2_f_NC))*vi2_ai2_Z_up

    res(:,FLZ) = (dlike(CL_f_NC) + dbarlike(CL_f_NC))*vi2_ai2_Z_down + &
         &       (ulike(CL_f_NC) + ubarlike(CL_f_NC))*vi2_ai2_Z_up

    res(:,F3Z) = (dlike(C3_f) - dbarlike(C3_f))*two_vi_ai_Z_down + &
         &       (ulike(C3_f) - ubarlike(C3_f))*two_vi_ai_Z_up
    res(:,F3Z) = res(:,F3Z)/xValues(grid)

    !--- deal with EM case -----------------------------------------
    res(:,F2EM) = (dlike(C2_f_NC) + dbarlike(C2_f_NC))*e2_down + &
         &        (ulike(C2_f_NC) + ubarlike(C2_f_NC))*e2_up
    
    res(:,FLEM) = (dlike(CL_f_NC) + dbarlike(CL_f_NC))*e2_down + &
         &        (ulike(CL_f_NC) + ubarlike(CL_f_NC))*e2_up

    ! for the W cases, it only makes sense to sum over an even number
    ! of light flavours; so save the actual number of flavours, switch
    ! the module-local cfh%nf_lcl variable to the nearest (lower) event number
    ! for our W calculations; switch back later
    nf_save = cfh%nf_lcl
    cfh%nf_lcl = (cfh%nf_lcl/2) * 2
    
    !--- deal with Wp case -----------------------------------------
    res(:,F2Wp) = (ulike(C2_f) + dbarlike(C2_f))*vi2_ai2_avg_W
    res(:,FLWp) = (ulike(CL_f) + dbarlike(CL_f))*vi2_ai2_avg_W
    res(:,F3Wp) = (ulike(C3_f) - dbarlike(C3_f))*two_vi_ai_avg_W
    res(:,F3Wp) = res(:,F3Wp)/xValues(grid)

    !--- deal with Wm case -----------------------------------------
    res(:,F2Wm) = (dlike(C2_f) + ubarlike(C2_f))*vi2_ai2_avg_W
    res(:,FLWm) = (dlike(CL_f) + ubarlike(CL_f))*vi2_ai2_avg_W
    res(:,F3Wm) = (dlike(C3_f) - ubarlike(C3_f))*two_vi_ai_avg_W
    res(:,F3Wm) = res(:,F3Wm)/xValues(grid)

    ! reset cfh%nf_lcl to the full (possibly odd) saved value
    cfh%nf_lcl = nf_save
    
    ! overall factor of two that we still haven't fully looked into as
    ! of 2015-02-24 [GPS TMP]
    res = res * two
    !! GPS+AK TMP: we have just included a factor of 2 but this plainly
    !! should not be present for the electromagnetic part; so here
    !! we eliminate it again...
    res(:,F2EM) = half * res(:,F2EM)
    res(:,FLEM) = half * res(:,FLEM)

  end function structure_function_general_full
  

  !----------------------------------------------------------------------
  function ulike(f) result(res)
    real(dp), intent(in) :: f(0:,ncompmin:)
    real(dp)             :: res(0:ubound(f,dim=1))
    res = sum(f(:, 2: cfh%nf_lcl: 2), dim=2)
  end function ulike
  !----------------------------------------------------------------------
  function dlike(f) result(res)
    real(dp), intent(in) :: f(0:,ncompmin:)
    real(dp)             :: res(0:ubound(f,dim=1))
    res = sum(f(:, 1: cfh%nf_lcl: 2), dim=2)
  end function dlike
  !----------------------------------------------------------------------
  function ubarlike(f) result(res)
    real(dp), intent(in) :: f(0:,ncompmin:)
    real(dp)             :: res(0:ubound(f,dim=1))
    res = sum(f(:,-2:-cfh%nf_lcl:-2), dim=2)
  end function ubarlike
  !----------------------------------------------------------------------
  function dbarlike(f) result(res)
    real(dp), intent(in) :: f(0:,ncompmin:)
    real(dp)             :: res(0:ubound(f,dim=1))
    res = sum(f(:,-1:-cfh%nf_lcl:-2), dim=2)
  end function dbarlike
  
  

  !----------------------------------------------------------------------
  real(dp) function alphasLocal(muR)
    real(dp), intent(in) :: muR
    real(dp) :: muR_lcl

    muR_lcl = max(muR,Qmin)
    ! we use alphas from the LHAPDF PDF
    ! alphasLocal = alphasPDF(muR_lcl)
    ! we use alphas from HOPPET
    alphasLocal = Value(coupling, muR_lcl)
  end function alphasLocal
  

  !----------------------------------------------------------------------
  ! F_LO
  ! calculate the leading order structure function at x, muF
  !
  function F_LO (y, Q) result(res)
    real(dp), intent(in)  :: y, Q
    real(dp) :: res(-6:6)
    real(dp) :: C1f(-6:6), C0P0f(-6:6)
    
    call EvalPdfTable_yQ(tables(1), y, Q, res)
    
  end function F_LO

  !----------------------------------------------------------------------
  ! F_NLO
  ! calculate the next-to-leading order structure function at x, muF
  !
  ! LRQ2 == ln muR^2/Q^2
  ! LFQ2 == ln muF^2/Q^2
  !
  function F_NLO (y, Q) result(res)
    real(dp), intent(in)  :: y, Q
    real(dp) :: res(-6:6), as2pi, LFQ2
    real(dp) :: C1f(-6:6), C0P0f(-6:6)

    as2pi = alphasLocal(Q * xmuR) / (twopi)

    ! C_NLO x f (x) in C1f(:) 
    call EvalPdfTable_yQ(tables(2), y, Q, C1f)
    res = C1f
    
    res = res * as2pi
    
  end function F_NLO


  !----------------------------------------------------------------------
  ! F_NNLO
  ! calculate the next-to-next-to-leading order structure function at x, muF
  !
  ! LRQ2 == ln muR^2/Q^2
  ! LFQ2 == ln muF^2/Q^2
  function F_NNLO (y, Q) result(res)
    real(dp), intent(in)  :: y, Q
    real(dp) :: res(-6:6), as2pi, LRQ2, LFQ2
    real(dp) :: C1f(-6:6), C0P0f(-6:6), C2f(-6:6), C0P0sqf(-6:6)
    real(dp) :: C0P1f(-6:6), C1P0f(-6:6)
    
    as2pi = alphasLocal(Q * xmuR) / (twopi)

    ! C_NNLO x f (x) in C2f(:,3) 
    call EvalPdfTable_yQ(tables(4), y, Q, C2f)
    res = C2f

    res = res * (as2pi)**2
    
  end function F_NNLO

  
  !----------------------------------------------------------------------
  ! F_N3LO
  ! calculate the next-to-next-to-next-to-leading order structure function at x, muF
  !
  ! LRQ2 == ln muR^2/Q^2
  ! LFQ2 == ln muF^2/Q^2
  function F_N3LO (y, Q) result(res)
    real(dp), intent(in)  :: y, Q
    real(dp) :: res(-6:6), as2pi, LRQ2, LFQ2
    real(dp) :: C1f(-6:6), C0P0f(-6:6), C3f(-6:6), C0P0sqf(-6:6), C2f(-6:6)
    real(dp) :: C0P1f(-6:6), C1P0f(-6:6), C0P0cbf(-6:6), C0P01f(-6:6), C0P10f(-6:6)
    real(dp) :: C1P0sqf(-6:6), C1P1f(-6:6), C2P0f(-6:6), C0P2f(-6:6)
    
    as2pi = alphasLocal(Q * xmuR) / (twopi)

    ! C_N3LO x f (x) in C2f(:,8) 
    call EvalPdfTable_yQ(tables(8), y, Q, C3f)
    res = C3f

    res = res * (as2pi)**3
    
  end function F_N3LO


  !----------------------------------------------------------------------
  ! dot product
  real(dp) function dot(p1,p2)
    real(dp), intent(in) :: p1(0:3), p2(0:3)
    dot = p1(0)*p2(0) - sum(p1(1:3)*p2(1:3))
  end function dot

  !----------------------------------------------------------------------
  ! set_scale_logs is only used for scale_choice = 0,1
  subroutine set_scale_logs(Q)
    real(dp), intent(in) :: Q

    log_muF2_over_Q2 = two * log(muF(Q) / Q)
    log_muR2_over_Q2 = two * log(muR(Q) / Q)
  end subroutine set_scale_logs
  
  !----------------------------------------------------------------------
  ! mu_R with scale_choice < 2
  real(dp) function muR(Q)
    real(dp), intent(in) :: Q
    muR = xmur * Q
  end function muR
  

!   !----------------------------------------------------------------------
!   ! mu_R1 as a function of Q1 and Q2
!   real(dp) function muR1(Q1, Q2, ptH)
!     real(dp), intent(in) :: Q1, Q2, ptH
!     muR1 = zero
!     if (scale_choice.le.1) then
!        ! if scale_choice = 0,1 then muR1(Q1,Q2) = muR(Q1)
!        muR1 = muR(Q1)
!     elseif (scale_choice.eq.2) then
!        ! else if scale_choice=2, use sqrt(Q1*Q2)
!        muR1 = xmur * sqrt(Q1 * Q2)
!     elseif (scale_choice.eq.3) then
!        ! else if scale_choice=3, use mixed scale
!        muR1 = xmur * mixed_scale(Q1, Q2, ptH)
!     endif
!   end function muR1
!   
!   !----------------------------------------------------------------------
!   ! mu_R2 as a function of Q1 and Q2
!   real(dp) function muR2(Q1, Q2, ptH)
!     real(dp), intent(in) :: Q1, Q2, ptH
!     muR2 = zero
!     if (scale_choice.le.1) then
!        ! if scale_choice = 0,1 then muR1(Q1,Q2) = muR(Q1)
!        muR2 = muR(Q2)
!     elseif (scale_choice.eq.2) then
!        ! else if scale_choice=2, use sqrt(Q1*Q2)
!        muR2 = xmur * sqrt(Q1 * Q2)
!     elseif (scale_choice.eq.3) then
!        ! else if scale_choice=3, use mixed scale
!        muR2 = xmur * mixed_scale(Q1, Q2, ptH)
!     endif
!   end function muR2

  !----------------------------------------------------------------------
  ! mu_F with scale_choice < 2
  real(dp) function muF(Q)
    real(dp), intent(in) :: Q
    muF = xmuf * Q
  end function muF
  
!   !----------------------------------------------------------------------
!   ! mu_F1 as a function of Q1 and Q2
!   real(dp) function muF1(Q1, Q2, ptH)
!     real(dp), intent(in) :: Q1, Q2, ptH
!     muF1 = zero
!     if (scale_choice.le.1) then
!        ! if scale_choice = 0,1 then muF1(Q1,Q2) = muF(Q1)
!        muF1 = muF(Q1)
!     elseif (scale_choice.eq.2) then
!        ! else if scale_choice=2, use sqrt(Q1*Q2)
!        muF1 = xmuf * sqrt(Q1 * Q2)
!     elseif (scale_choice.eq.3) then
!        ! else if scale_choice=3, use mixed scale
!        muF1 = xmuf * mixed_scale(Q1, Q2, ptH)
!     else
!        call wae_error('muF1(Q)', 'illegal value for scale_choice', intval = scale_choice)
!     endif
!   end function muF1
! 
!   !----------------------------------------------------------------------
!   ! mu_F2 as a function of Q1 and Q2
!   real(dp) function muF2(Q1, Q2, ptH)
!     real(dp), intent(in) :: Q1, Q2, ptH
!     muF2 = zero
!     if (scale_choice.le.1) then
!        ! if scale_choice = 0,1 then muF1(Q1,Q2) = muF(Q1)
!        muF2 = muF(Q2)
!     elseif (scale_choice.eq.2) then
!        ! else if scale_choice=2, use sqrt(Q1*Q2)
!        muF2 = xmuf * sqrt(Q1 * Q2)
!     elseif (scale_choice.eq.3) then
!        ! else if scale_choice=3, use mixed scale
!        muF2 = xmuf * mixed_scale(Q1, Q2, ptH)
!     else
!        call wae_error('muF2(Q)', 'illegal value for scale_choice', intval = scale_choice)
!     endif
!   end function muF2

!   !----------------------------------------------------------------------
!   ! Defines which scale to use for scale_choice = 3,
!   ! which can be any function of Q1, Q2 and ptH
!   real(dp) function mixed_scale(Q1,Q2,ptH)
!     real(dp), intent(in) :: Q1, Q2, ptH
!     mixed_scale = ((mh*0.5d0)**4d0+(mh*ptH*0.5d0)**2d0)**(0.25d0)
!   end function mixed_scale

end module hoppet_structure_functions
