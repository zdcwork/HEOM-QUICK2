subroutine check_inum_vec(inp, nvec, ndim, ipost, iout)
implicit none
!
! purpose: check whether input number exists in the vector
!          output its first position (npost)
!          iout = 0 if the input exists in the vector
!
integer, intent(in)  :: inp, ndim
integer, intent(in)  :: nvec(*)
integer, intent(out) :: ipost, iout
integer              :: ni
!
ipost = 0
iout  = -1
!
do ni=1,ndim
   if (inp .eq. nvec(ni)) then
       ipost = ni
       iout  = 0
       exit 
   end if
end do
!
return
end subroutine check_inum_vec
!
subroutine comp_inum_vec12(nvec1, ndim1, nvec2, ndim2, nvec_tmp1, nvec_tmp2, imode, ipost)
implicit none
!
! purpose: compare two integer vectors (not necessarily ordered)
!
!   input: imode=1: find first element in nvec1 that is not in nvec2
!          imode=2: find first element in nvec2 that is not in nvec1
!               otherwise stop
!
!  output: ipost: the position of the found element, 0 if not found
! 
integer, intent(in) :: ndim1, ndim2, imode
integer, intent(in) :: nvec1(*), nvec2(*)
integer, intent(inout) :: nvec_tmp1(*), nvec_tmp2(*)
integer, intent(out) :: ipost
integer              :: ni, nj, nk
!
ipost = 0
!
if (imode .eq. 1) then
    nvec_tmp2(1:ndim2) = nvec2(1:ndim2)
    do ni=1,ndim1
       nk = 1
       loopj1: do nj=1,ndim2
          if (nvec1(ni) .eq. nvec_tmp2(nj)) then
              nvec_tmp2(nj) = -1
              exit loopj1
          else
              nk = nk + 1
          end if
       end do loopj1
       if (nk .gt. ndim2) then
           ipost = ni
           return
       end if
    end do 
!
else if (imode .eq. 2) then
    nvec_tmp1(1:ndim1) = nvec1(1:ndim1)
    do ni=1,ndim2
       nk = 1
       loopj2: do nj=1,ndim1
          if (nvec2(ni) .eq. nvec_tmp1(nj)) then
              nvec_tmp1(nj) = -1
              exit loopj2
          else
              nk = nk + 1
          end if
       end do loopj2
       if (nk .gt. ndim1) then
           ipost = ni
           return
       end if
    end do 
!
else
    write(6,*)
    write(6,*)'comp_inum_vec12: error! unknown imode ', imode
    write(6,*)' nvec1:  ndim1 =  ', ndim1
    write(6,*)(nvec1(ni), ni=1,ndim1)
    write(6,*)' nvec2:  ndim2 =  ', ndim2
    write(6,*)(nvec2(ni), ni=1,ndim2)
    stop
end if
! 
return
end subroutine comp_inum_vec12
