!--------------------------------------------------------------------
! contains the random unmber generator; for the time being keep it
! simple and just use the built-in one -- then replace it with
! something better shortly
module random
  use types, only : dp
  implicit none

  !-- the following combination causes problems with pgf90
  !   so eliminate it, and shield the dp stuff
  !private
  !public :: ran,rangen
  private :: dp

contains
  function ran()
    real(dp) :: ran
    real(dp) :: rann(1)
    !call random_number(ran)
    !call rangen(1,rann)
    call rm48(rann,1)
    ran = rann(1)
  end function ran
 
  SUBROUTINE RANGEN(N,R)
    !IMPLICIT NONE
    !---RANDOM NUMBER GENERATOR
    !   USES METHOD OF l'Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
    !   RETURNS A VECTOR OF N RANDOM VALUES
    !   IF (N.EQ.0) THE FIRST TWO VALUES IN R SET THE SEEDS
    !   IF (N.LT.0) PRINT THE CURRENT VALUES OF THE SEEDS
    real(dp) r(*)
      integer n,iseed(3)
      if (n.gt.0) then
         call rm48(r,n)
      else if (n.lt.0) then
         call rm48ut(iseed(1),iseed(2),iseed(3))
         write (*,'(i10,a,i10,i11,i11)') -n-1,', iseed=',&
              &        iseed(1),iseed(2),iseed(3)
      else ! n=0
        if(nint(r(1)) .eq. 0) then
           !-- retrieve seed --
           call rm48ut(iseed(1),iseed(2),iseed(3))
           r(1) = iseed(1)
           r(2) = iseed(2)
           r(3) = iseed(3)
        else
           !-- set seed -------
           iseed(1)=nint(r(1))
           iseed(2)=nint(r(2))
           iseed(3)=nint(r(3))
           call rm48in(iseed(1),iseed(2),iseed(3))
        end if
      end if
  END SUBROUTINE RANGEN
 
end module random


