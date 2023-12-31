subroutine genpermute(nlen, vecinp, vecout, fact)
implicit none
!
! generate a new permutation, vecout, from an old permutation, vecinp,
! in lexicographic order
!
integer, intent(in)  :: nlen, vecinp(*)
integer, intent(out) :: vecout(*), fact
!
integer, parameter   :: mone = -1
!
integer              :: ni, nj, nk, nl
integer              :: ntmp1, ntmp2, ntmp3
!
if (nlen .le. 1) then
    write(6,*)
    write(6,*)'genpermute: error! nlen = ', nlen
    stop
end if
!
vecout(1:nlen) = vecinp(1:nlen)
fact = 1
!
! find nj
! 
nj = nlen - 1
do while (nj .ge. 1 .and. vecout(nj) .ge. vecout(nj+1))
   nj = nj - 1
end do
!
if (nj .eq. 0) then
    write(6,*)
    write(6,*)'genpermute: error! vecinp is already the last permutation '
    write(6,*) (vecinp(ni), ni=1,nlen)
    write(6,*) nj
    stop
end if
!
! find nl
! 
nl = nlen
do while (vecout(nj) .ge. vecout(nl))
   nl = nl - 1
end do
!
! exchange vecout(nj) and vecout(nl)
!
if (nj .ne. nl) then
    ntmp1      = vecout(nj)
    vecout(nj) = vecout(nl)
    vecout(nl) = ntmp1
    fact       = fact * mone
end if
!
! revert sub-vector vecout(nj+1,nlen)
!
ntmp1 = nlen - nj
ntmp2 = ntmp1 / 2
ntmp3 = ntmp1 * (ntmp1 - 1) / 2
fact  = fact * mone**ntmp3
!
if (ntmp1 .gt. 1) then
    do ni=1,ntmp2
       ntmp1 = vecout(nj+ni)
       vecout(nj+ni) = vecout(nlen-ni+1)
       vecout(nlen-ni+1) = ntmp1
    end do
end if
!
return
end subroutine genpermute
