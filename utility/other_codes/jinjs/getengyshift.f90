subroutine getengyshift(tt, lead, ispin, eshift)
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in)  :: lead, ispin     ! 1 for left, 2 for right
real*8,  intent(in)  :: tt
real*8,  intent(out) :: eshift
!
if (lead .ne. 1 .and. lead .ne. 2) then
  write(6,*)
  write(6,*)' error! unknown lead in <getengyshift>! ', lead
  stop
end if
!
if (lead .eq. 1 .and. tt .gt. toffL) then
  eshift = 0.d0
  return
else if (lead .eq. 2 .and. tt .gt. toffR) then
  eshift = 0.d0
  return
end if
!
if (fieldtype .eq. 0) then
  if (lead .eq. 1) then
    if (.not. leadspin) then
      eshift = ampL * (1.d0 - dexp(-tt / aL))
    else
      eshift = engyshift(1,ispin) * (1.d0 - dexp(-tt / aL))
    end if
  else
    eshift = ampR * (1.d0 - dexp(-tt / aR))
  end if
else if (fieldtype .eq. 1) then
  if (lead .eq. 1) then
    eshift = ampL * (1.d0 - dcos(wL * tt))
  else
    eshift = ampR * (1.d0 - dcos(wR * tt))
  end if
else if (fieldtype .eq. 2) then 
  if (lead .eq. 1) then
    if (tt .lt. tuL) then
      eshift = ampL * (1.d0 - dexp(-tt / auL))
    else if (tt .ge. tuL .and. tt .lt. tdL) then
      eshift = ampL * (1.d0 - dexp(-tuL / auL))
    else 
      eshift = ampL * (1.d0 - dexp(-tuL / auL)) * dexp(-(tt - tdL) / adL)
    end if
  else
    if (tt .lt. tuR) then
      eshift = ampR * (1.d0 - dexp(-tt / auR))
    else if (tt .ge. tuR .and. tt .lt. tdR) then
      eshift = ampR * (1.d0 - dexp(-tuR / auR))
    else 
      eshift = ampR * (1.d0 - dexp(-tuR / auR)) * dexp(-(tt - tdR) / adR)
    end if
  end if
else if (fieldtype .eq. 3) then
  if (lead .eq. 1) then
    if (tt .le. tcl) then
      eshift = ampL * dexp( -((tt - tcl)/taul)**2 )
    else
      eshift = ampL * dexp( -( (tt - tcl)/(taul/2.d0) )**2 )
!      eshift = ampL * dexp( -( (tt - tcl)/(taul/1.d0) )**2 )
    end if
  else 
    if (tt .le. tcr) then
      eshift = ampR * dexp( -((tt - tcr)/taur)**2 )
    else
      eshift = ampR * dexp( -( (tt - tcr)/(taur/2.d0) )**2 )
!      eshift = ampR * dexp( -( (tt - tcr)/(taur/1.d0) )**2 )
    end if
  end if
else if (fieldtype .eq. 4) then
  if (lead .eq. 1) then
    if (tt .lt. tonL) then
      eshift = 0.d0
    else if (tt .ge. tonL .and. tt .lt. tholdL) then
      eshift = (tt - tonL) / (tholdL - tonL) * ampL
    else
      eshift = ampL
    end if
  else
    if (tt .lt. tonR) then
      eshift = 0.d0
    else if (tt .ge. tonR .and. tt .lt. tholdR) then
      eshift = (tt - tonR) / (tholdR - tonR) * ampR
    else
      eshift = ampR
    end if
  end if
else
  write(6,*)
  write(6,*)' error! unknown fieldtype in <getengyshift>! ', fieldtype
  stop
end if
!
end subroutine getengyshift
