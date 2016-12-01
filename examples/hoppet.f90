
module hoppet
  use types
  use qed_coupling_module
  implicit none
  type(qed_coupling) :: alpha_qed
  logical :: initialized = .false.

contains
  subroutine initialize(mcharm, mbottom, mtop)
    real(dp), intent(in) :: mcharm, mbottom, mtop
    call InitQEDCoupling(alpha_qed, mcharm, mbottom, mtop)
    initialized = .true.
  end subroutine

  real(dp) function alphaqed(Q) result(res)
    real(dp), intent(in) :: Q
    if (.not.initialized) then
      call exit()
    end if
    res = Value(alpha_qed, Q)
  end function alphaqed

end module hoppet
