
module hoppet
  use types
  use qed_coupling_module
  use qcd_coupling
  use hoppet_structure_functions
  implicit none
  type(qed_coupling) :: alpha_qed
  logical :: initialized = .false.
  real(dp) :: quark_masses(1:6) = -1d0

contains
  subroutine initialize(mcharm, mbottom, mtop)
    real(dp), intent(in) :: mcharm, mbottom, mtop
    integer :: iflv
    call InitQEDCoupling(alpha_qed, mcharm, mbottom, mtop)
    call InitPDFsetByName("PDF4LHC15_nnlo_100")
    call InitPDFm(1,0)
    call hoppet_start()
    call hoppet_read_PDF()
    call hoppet_set_structure_functions_up_to(3)
    do iflv = 4, 6
      quark_masses(iflv) = QuarkMass(coupling, iflv)
    end do
    initialized = .true.
  end subroutine

  subroutine checkinit()
    if (.not.initialized) then
      print *, "Error HOPPET not initialized"
      call exit()
    end if
  end subroutine

  real(dp) function alphaqed(Q) result(res)
    real(dp), intent(in) :: Q
    call checkinit()
    res = Value(alpha_qed, Q)
  end function alphaqed

  real(dp) function f2(x,Q) result(res)
    real(dp), intent(in) :: x, Q
    real(dp) :: Fi(-6:6), y_ln1ox
    call checkinit()
    y_ln1ox = log(1.0d0/x)
    Fi =      F_LO  (y_ln1ox, Q)
    Fi = Fi + F_NLO (y_ln1ox, Q)
    Fi = Fi + F_NNLO(y_ln1ox, Q)
    res = Fi(F2EM)
  end function f2

  real(dp) function fl(x,Q) result(res)
    real(dp), intent(in) :: x, Q
    real(dp) :: Fi(-6:6), y_ln1ox
    call checkinit()
    y_ln1ox = log(1.0d0/x)
    Fi =      F_LO  (y_ln1ox, Q)
    Fi = Fi + F_NLO (y_ln1ox, Q)
    Fi = Fi + F_NNLO(y_ln1ox, Q)
    res = Fi(FLEM)
  end function fl

  real(dp) function masses(fl) result(res)
    integer, intent(in) :: fl
    res = quark_masses(fl)
  end function masses

end module hoppet
