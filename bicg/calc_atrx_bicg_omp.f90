subroutine calc_atrx_bicg_omp(dim0, rcgtmp, rcgrhs)
!
! purpose : calculate the RHS of EOM for rho
!           and replace the first equation with the normalization 
!           condition for system density matrix.
!    
! input/output : rcgtmp, rcgrhs (refreshed)
!
! note :
!       Here, the Hermicity of rcgtmp(*,*,1) has been considered by calling
!        the subroutine <modifylmat>, where the Liouville coefficient matrix
!        has been rearranged.
!
!       To ignore a certain unknown variable from the VECTOR rcgtmp, just
!        assign zero value to the corresponding position, and also put zero 
!        at the same position of the VECTOR rcgrhs.        
!       However, to guarantee a correct Liouvill coefficient matrix can be 
!        achieved, such as that constructed in <checkatrx.f90>, corresponding
!        Liouville matrices need to be modified (remove lines and columns from
!        A^t). For instance, see subroutine <modifylmat>.
!
use omp_lib
use auxmod
use auxmod_omp
use matmod
use matmod_omp
use tmpmatmod
use tmpmatmod_omp
use tmplioumat_omp
use bicgmod_omp
implicit none
include '../include/sizes'
include '../include/common'
!
integer,    intent(in)    :: dim0
complex*16, intent(inout) :: rcgtmp(dim0,dim0,*), rcgrhs(dim0,dim0,*)
integer    :: itier, iorbs, iorbs2, ispin, ialf, isgn, imats, jsgn, itype
integer    :: ifactfront, ifactrear
integer    :: index_nk(MAXTIER, maxprocs)
integer*8  :: lni, lnj, lnl, lnm
integer    :: ni, nj, nk, nl, nm, nn
integer    :: ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, dim1
integer    :: i, j, k
real*8     :: cpu1, cpu2, dtmp1, dtmp2
real*8     :: dmaxhi, dbr, dbi, ddr, ddi
logical    :: lhsi
complex*16 :: cgamma_nk
integer    :: itimes
data itimes /0/
!
integer    :: icore
!
if (itimes .eq. 0) then
   write(6,*)'calc_atrx_bicg_omp: first entry'
end if
!
itimes = itimes + 1
!
if (dim0 .ne. nrho) then
  write(6,*)'calc_atrx_bicg_omp: error! dim0 != nrho ', dim0, nrho
  stop
end if
dim1 = nrho**2
!
if (lhf) then
   call calchs(rcgtmp(1,1,1))
end if
if (igroundsteady .eq. 0) then
   cmtmp3 = hs
else if (igroundsteady .eq. 1) then
   call calcdhs(1.d20, rcgtmp(1,1,1))
   cmtmp3 = hs + dhs
else
  write(6,*)
  write(6,*)' error! unknown igroundsteady in calc_atrx_bicg_omp ', igroundsteady
  stop
end if
!
dmtmp1 =  dble( -eye * cmtmp3 )
dmtmp2 = dimag( -eye * cmtmp3 )
call dmaxmat(nrho, nrho, dmtmp2, nrho, dmaxhi)
!
call h2lmln(nrho, dmtmp1, nrho, lhsrnl, dim1)
call h2lmrn(nrho, dmtmp1, nrho, lhsrnr, dim1)
if (dmaxhi .ge. dnano) then
  lhsi = .true.
  call h2lmln(nrho, dmtmp2, nrho, lhsinl, dim1)
  call h2lmrn(nrho, dmtmp2, nrho, lhsinr, dim1)
else
  lhsi   = .false.
  lhsinl = 0.d0
  lhsinr = 0.d0
end if
!
!lnl = 1
!
!######################################################################################
!                                                                                     !
!                   Follow the sequence settled in indextable                         !
!                                                                                     !
!######################################################################################
!
do ni=1,nrho
  rcgtmp(ni,ni,1) = dcmplx(dble(rcgtmp(ni,ni,1)), 0.d0)
end do
do ni=1,nrho
  do nj=ni+1,nrho
    rcgtmp(ni,nj,1) = czero
  end do
end do
rcgrhs(1:nrho, 1:nrho, 1:nunk) = czero
!
! itier=1, nballs=0,  no balls, i.e., N = 0 in the derivation
!
itier = 1
dmtmp1(1:nrho, 1:nrho) =  dble(rcgtmp(1:nrho, 1:nrho, 1))
dmtmp2(1:nrho, 1:nrho) = dimag(rcgtmp(1:nrho, 1:nrho, 1))
call h2lvec(nrho, dmtmp1, nrho, dlvec1) 
call h2lvec(nrho, dmtmp2, nrho, dlvec2) 
!
dlmat1 = lhsrnl - lhsrnr
call modifylmat(nrho, dlmat1, dim1, 2)
dlmat3 = dlmat1
call modifylmat(nrho, dlmat1, dim1, 0)
call modifylmat(nrho, dlmat3, dim1, 1)
call dgemv('t', dim1, dim1, 1.d0, dlmat1, dim1, dlvec1, 1, 0.d0, dlvec3, 1)
call dgemv('t', dim1, dim1, 1.d0, dlmat3, dim1, dlvec2, 1, 0.d0, dlvec4, 1)
if (lhsi) then
  dlmat2 = lhsinl - lhsinr
  call modifylmat(nrho, dlmat2, dim1, 2)
  dlmat4 = dlmat2
  call modifylmat(nrho, dlmat2, dim1, 0)
  call modifylmat(nrho, dlmat4, dim1, 1)
  call dgemv('t', dim1, dim1, 1.d0, dlmat2, dim1, dlvec2, 1, 1.d0, dlvec3, 1)
  call dgemv('t', dim1, dim1,-1.d0, dlmat4, dim1, dlvec1, 1, 1.d0, dlvec4, 1)
end if
call l2hvec(nrho, dmtmp5, nrho, dlvec3)
call l2hvec(nrho, dmtmp6, nrho, dlvec4)
rcgrhs(1:nrho, 1:nrho, 1) = rcgrhs(1:nrho, 1:nrho, 1) +                            &
                            dcmplx(dmtmp5(1:nrho, 1:nrho), dmtmp6(1:nrho, 1:nrho))
!
!==========================================
!$omp parallel default(shared) private(ni,nj,isgn,jsgn,ialf,iorbs,ispin,imats,dtmp1,dtmp2,lni,lnj,nn,icore)
icore = omp_get_thread_num() + 1
dlvec1_omp(1:dim1,icore) = dlvec1(1:dim1)
dlvec2_omp(1:dim1,icore) = dlvec2(1:dim1)
!
!$omp do schedule(dynamic)
do lnj=ifirst(2), ilast(2)
  nj = indextable(lnj)
  if (mpm(nj) .ne. 1) then
    write(6,*)
    write(6,*)' error! wrong mpm for itier=2 ', lni, nj, mpm(nj)
    stop
  end if
  isgn   = 1
  jsgn   = 2
  ialf   = ilead(nj)
  iorbs  = morbs(nj)
  ispin  = mspin(nj)
  imats  = mcor(nj)
!  
  if (lscale) then
     dtmp1 = dbsqrt(iorbs,ispin,imats,ialf,isgn)
     dtmp2 = dbsqrt(iorbs,ispin,imats,ialf,jsgn)
  else
     dtmp1 = 1.d0
     dtmp2 = 1.d0
  end if
!
  call look4drawer(isgn, ialf, iorbs, ispin, imats, nn)
  lni = nfirst(2) - 1 + nn
  dlvec7_omp(1:dim1,icore) = 0.d0   ! Yr
  dlvec8_omp(1:dim1,icore) = 0.d0   ! Yi
! -i * c * rho
  call lcpams_omp (iorbs, ispin, jsgn, 1, 1, icore)
  call modlams_omp(iorbs, ispin, jsgn, 1, 1, 2, icore)
  call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                    lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec3_omp(1,icore))
  call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                    lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec4_omp(1,icore))
  dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) - dlvec3_omp(1:dim1,icore) * dtmp1
  dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) + dlvec4_omp(1:dim1,icore) * dtmp1
! i * rho * c
  call lcpams_omp (iorbs, ispin, jsgn, 2, 1, icore)
  call modlams_omp(iorbs, ispin, jsgn, 2, 1, 2, icore)
  call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                    lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec3_omp(1,icore))
  call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                    lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec4_omp(1,icore))  
  dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + dlvec3_omp(1:dim1,icore) * dtmp1
  dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) - dlvec4_omp(1:dim1,icore) * dtmp1
! -i * c^t * rho^dag
  call lcpams_omp (iorbs, ispin, isgn, 1, 2, icore)
  call modlams_omp(iorbs, ispin, isgn, 1, 2, 2, icore)
  call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                    lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec3_omp(1,icore))
  call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                    lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec4_omp(1,icore))  
  dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) - dlvec3_omp(1:dim1,icore) * dtmp2
  dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) - dlvec4_omp(1:dim1,icore) * dtmp2
! i * rho^dag * c^t
  call lcpams_omp (iorbs, ispin, isgn, 2, 2, icore)
  call modlams_omp(iorbs, ispin, isgn, 2, 2, 2, icore)
  call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                    lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec3_omp(1,icore))
  call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                    lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec4_omp(1,icore))  
  dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + dlvec3_omp(1:dim1,icore) * dtmp2
  dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) + dlvec4_omp(1:dim1,icore) * dtmp2
!
  call l2hvec(nrho, dmtmp5_omp(1,1,icore), nrho, dlvec7_omp(1,icore))
  call l2hvec(nrho, dmtmp6_omp(1,1,icore), nrho, dlvec8_omp(1,icore))
  rcgrhs(1:nrho, 1:nrho, lni) = rcgrhs(1:nrho, 1:nrho, lni) +                         &
                     dcmplx(dmtmp5_omp(1:nrho, 1:nrho, icore), dmtmp6_omp(1:nrho, 1:nrho, icore))
end do
!$omp end do
!$omp end parallel
!============================================
!rcgrhs_omp(1:nrho,1:nrho,1:nunk,1:nthreads_omp) = czero
!
! Higher tiers, # of balls = itier - 1. Here, itier >= 2, i.e., N >= 1.
!
do itier=2,ntier
  call cpu_time(cpu1)

!$omp parallel default(shared) private(lni,lnj,lnl,cgamma_nk,ni,nj,itype,ifactfront,ifactrear,nn,isgn,ialf,iorbs) &
!$omp private(ispin,imats,iorbs2,dbr,dbi,ddr,ddi,ntmp5,jsgn,dtmp1,dtmp2,icore)
  icore = omp_get_thread_num() + 1
!$omp do schedule(dynamic)
  do lni=nfirst(itier), nlast(itier)
!
    cgamma_nk = czero
!
! Get the multi-component index {n_k}, where k goes through 1, ..., nvar.   
!
    do ni=1,itier-1 
     index_nk(ni,icore) = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
     cgamma_nk    = cgamma_nk + cgamma(index_nk(ni,icore))
     if (igroundsteady .ne. 0) then
        cgamma_nk = cgamma_nk + eye * dipm(index_nk(ni,icore)) *                    &
                    eleadinfty(ilead(index_nk(ni,icore)), mspin(index_nk(ni,icore)))
     end if
    end do
!
! same itier contribution (N)
!
    dmtmp1_omp(1:nrho, 1:nrho, icore) =  dble(rcgtmp(1:nrho, 1:nrho, lni))
    dmtmp2_omp(1:nrho, 1:nrho, icore) = dimag(rcgtmp(1:nrho, 1:nrho, lni))
    call h2lvec(nrho, dmtmp1_omp(1,1,icore), nrho, dlvec1_omp(1,icore)) 
    call h2lvec(nrho, dmtmp2_omp(1,1,icore), nrho, dlvec2_omp(1,icore)) 
!
    cmtmp1_omp(1:nrho,1:nrho,icore) = -eye * cmtmp3(1:nrho,1:nrho)
    if (lhsi) then
      call lioumatmat(1, 0, 1, 1, nrho, cmtmp1_omp(1,1,icore), nrho, rcgtmp(1,1,lni), nrho,  czero, cmtmp2_omp(1,1,icore), nrho)
      cmtmp1_omp(1:nrho,1:nrho,icore) = -cmtmp1_omp(1:nrho,1:nrho,icore)
      call lioumatmat(2, 0, 1, 1, nrho, cmtmp1_omp(1,1,icore), nrho, rcgtmp(1,1,lni), nrho, cunity, cmtmp2_omp(1,1,icore), nrho)
    else
      call lioumatmat(1, 0, 0, 1, nrho, cmtmp1_omp(1,1,icore), nrho, rcgtmp(1,1,lni), nrho,  czero, cmtmp2_omp(1,1,icore), nrho)
      cmtmp1_omp(1:nrho,1:nrho,icore) = -cmtmp1_omp(1:nrho,1:nrho,icore)
      call lioumatmat(2, 0, 0, 1, nrho, cmtmp1_omp(1,1,icore), nrho, rcgtmp(1,1,lni), nrho, cunity, cmtmp2_omp(1,1,icore), nrho)
    end if
    rcgrhs(1:nrho, 1:nrho, lni) = rcgrhs(1:nrho, 1:nrho, lni) + cmtmp2_omp(1:nrho, 1:nrho, icore) + &
                                  dconjg(cgamma_nk) * rcgtmp(1:nrho, 1:nrho, lni)
!
    lnl = index_coef_ref(lni)
!
! nearest lower tier contribution (N - 1)
!
    do ni=1,itier-1
      lnj        = indexcoef(lnl)
      itype      = itypecoef(lnl)
      ifactfront = ifactfcoef(lnl)
      ifactrear  = ifactrcoef(lnl)
      lnl        = lnl + 1
!----
!cycle ! -
!----
      nn    = index_nk(ni,icore)
      isgn  = mpm(nn)
      ialf  = ilead(nn) 
      iorbs = morbs(nn)
      ispin = mspin(nn)
      imats = mcor(nn)
!
      dlvec7_omp(1:dim1,icore) = 0.d0       ! Yr
      dlvec8_omp(1:dim1,icore) = 0.d0       ! Yi
!
      if (.not. offcor) then
        iorbs2 = iorbs
        dbr    =  dble( -eye * cb(iorbs,ispin,imats,ialf,isgn) * dble(ifactfront) )
        dbi    = dimag( -eye * cb(iorbs,ispin,imats,ialf,isgn) * dble(ifactfront) )
        ddr    =  dble(  eye * cd(iorbs,ispin,imats,ialf,isgn) * dble(ifactrear)  )
        ddi    = dimag(  eye * cd(iorbs,ispin,imats,ialf,isgn) * dble(ifactrear)  )
      end if
!
      if (.not. offcor) then
        if (itype .eq. 0) then    ! rcgtmp
          call lcpams_omp (iorbs2, ispin, isgn, 1, 1, icore)
          if (itier .eq. 2) call modlams_omp(iorbs2, ispin, isgn, 1, 1, 0, icore)
          call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                            lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec3_omp(1,icore))
          dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + dbr * dlvec3_omp(1:dim1,icore)
          call lcpams_omp (iorbs2, ispin, isgn, 1, 1, icore)
          if (itier .eq. 2) call modlams_omp(iorbs2, ispin, isgn, 1, 1, 1, icore)
          call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                            lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec4_omp(1,icore))
          dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) + dbr * dlvec4_omp(1:dim1,icore)
          if (dabs(dbi) .ge. dnano) then
            if (itier .ne. 2) then
              dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + dbi * dlvec4_omp(1:dim1,icore)
              dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) - dbi * dlvec3_omp(1:dim1,icore) 
            else
              call lcpams_omp (iorbs2, ispin, isgn, 1, 1, icore)
              call modlams_omp(iorbs2, ispin, isgn, 1, 1, 0, icore)
              call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                                lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec4_omp(1,icore))
              dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + dbi * dlvec4_omp(1:dim1,icore)
              call lcpams_omp (iorbs2, ispin, isgn, 1, 1, icore)
              call modlams_omp(iorbs2, ispin, isgn, 1, 1, 1, icore)
              call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                                lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec3_omp(1,icore))
              dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) - dbi * dlvec3_omp(1:dim1,icore)
            end if
          end if
!
          call lcpams_omp (iorbs2, ispin, isgn, 2, 1, icore)
          if (itier .eq. 2) call modlams_omp(iorbs2, ispin, isgn, 2, 1, 0, icore)
          call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                            lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec3_omp(1,icore))
          dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + ddr * dlvec3_omp(1:dim1,icore)
          call lcpams_omp (iorbs2, ispin, isgn, 2, 1, icore)
          if (itier .eq. 2) call modlams_omp(iorbs2, ispin, isgn, 2, 1, 1, icore)
          call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                            lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec4_omp(1,icore))
          dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) + ddr * dlvec4_omp(1:dim1,icore)
          if (dabs(ddi) .ge. dnano) then
            if (itier .ne. 2) then
              dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + ddi * dlvec4_omp(1:dim1,icore)
              dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) - ddi * dlvec3_omp(1:dim1,icore) 
            else
              call lcpams_omp (iorbs2, ispin, isgn, 2, 1, icore)
              call modlams_omp(iorbs2, ispin, isgn, 2, 1, 0, icore)
              call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                                lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec4_omp(1,icore))
              dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + ddi * dlvec4_omp(1:dim1,icore)
              call lcpams_omp (iorbs2, ispin, isgn, 2, 1, icore)
              call modlams_omp(iorbs2, ispin, isgn, 2, 1, 1, icore)
              call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                                lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec3_omp(1,icore))
              dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) - ddi * dlvec3_omp(1:dim1,icore)
            end if
          end if
        else                      ! rcgtmp^dag
          call lcpams_omp (iorbs2, ispin, isgn, 1, 2, icore)
          call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                            lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec3_omp(1,icore))
          dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + dbr * dlvec3_omp(1:dim1,icore)
          call lcpams_omp (iorbs2, ispin, isgn, 1, 2, icore)
          call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                            lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec4_omp(1,icore))
          dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) - dbr * dlvec4_omp(1:dim1,icore)
          if (dabs(dbi) .ge. dnano) then
            dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + dbi * dlvec4_omp(1:dim1,icore)
            dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) + dbi * dlvec3_omp(1:dim1,icore) 
          end if
!
          call lcpams_omp (iorbs2, ispin, isgn, 2, 2, icore)
          call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                            lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec3_omp(1,icore))
          dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + ddr * dlvec3_omp(1:dim1,icore)
          call lcpams_omp (iorbs2, ispin, isgn, 2, 2, icore)
          call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                            lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec4_omp(1,icore))
          dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) - ddr * dlvec4_omp(1:dim1,icore)
          if (dabs(ddi) .ge. dnano) then
            dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + ddi * dlvec4_omp(1:dim1,icore)
            dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) + ddi * dlvec3_omp(1:dim1,icore) 
          end if
        end if
!
      else ! offcor
        do iorbs2=1,norbs
          ntmp5 = (iorbs2 - 1) * norbs + iorbs
          dbr   =  dble( -eye * cb(ntmp5,ispin,imats,ialf,isgn) * dble(ifactfront) )
          dbi   = dimag( -eye * cb(ntmp5,ispin,imats,ialf,isgn) * dble(ifactfront) )
          ddr   =  dble(  eye * cd(ntmp5,ispin,imats,ialf,isgn) * dble(ifactrear)  )
          ddi   = dimag(  eye * cd(ntmp5,ispin,imats,ialf,isgn) * dble(ifactrear)  )
          if (itype .eq. 0) then    ! rcgtmp
            call lcpams_omp (iorbs2, ispin, isgn, 1, 1, icore)
            if (itier .eq. 2) call modlams_omp(iorbs2, ispin, isgn, 1, 1, 0, icore)
            call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                              lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec3_omp(1,icore))
            dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + dbr * dlvec3_omp(1:dim1,icore)
            call lcpams_omp (iorbs2, ispin, isgn, 1, 1, icore)
            if (itier .eq. 2) call modlams_omp(iorbs2, ispin, isgn, 1, 1, 1, icore)
            call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                              lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec4_omp(1,icore))
            dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) + dbr * dlvec4_omp(1:dim1,icore)
            if (dabs(dbi) .ge. dnano) then
              if (itier .ne. 2) then
                dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + dbi * dlvec4_omp(1:dim1,icore)
                dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) - dbi * dlvec3_omp(1:dim1,icore)
              else
                call lcpams_omp (iorbs2, ispin, isgn, 1, 1, icore)
                call modlams_omp(iorbs2, ispin, isgn, 1, 1, 0, icore)
                call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                                  lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec4_omp(1,icore))
                dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + dbi * dlvec4_omp(1:dim1,icore)
                call lcpams_omp (iorbs2, ispin, isgn, 1, 1, icore)
                call modlams_omp(iorbs2, ispin, isgn, 1, 1, 1, icore)
                call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                                  lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec3_omp(1,icore))
                dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) - dbi * dlvec3_omp(1:dim1,icore)
              end if
            end if
!
            call lcpams_omp (iorbs2, ispin, isgn, 2, 1, icore)
            if (itier .eq. 2) call modlams_omp(iorbs2, ispin, isgn, 2, 1, 0, icore)
            call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                              lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec3_omp(1,icore))
            dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + ddr * dlvec3_omp(1:dim1,icore)
            call lcpams_omp (iorbs2, ispin, isgn, 2, 1, icore)
            if (itier .eq. 2) call modlams_omp(iorbs2, ispin, isgn, 2, 1, 1, icore)
            call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                              lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec4_omp(1,icore))
            dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) + ddr * dlvec4_omp(1:dim1,icore)
            if (dabs(ddi) .ge. dnano) then
              if (itier .ne. 2) then
                dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + ddi * dlvec4_omp(1:dim1,icore)
                dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) - ddi * dlvec3_omp(1:dim1,icore) 
              else
                call lcpams_omp (iorbs2, ispin, isgn, 2, 1, icore)
                call modlams_omp(iorbs2, ispin, isgn, 2, 1, 0, icore)
                call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                                  lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec4_omp(1,icore))
                dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + ddi * dlvec4_omp(1:dim1,icore)
                call lcpams_omp (iorbs2, ispin, isgn, 2, 1, icore)
                call modlams_omp(iorbs2, ispin, isgn, 2, 1, 1, icore)
                call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                                  lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec3_omp(1,icore))
                dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) - ddi * dlvec3_omp(1:dim1,icore)
              end if
            end if
          else                      ! rcgtmp^dag
            call lcpams_omp (iorbs2, ispin, isgn, 1, 2, icore)
            call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                              lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec3_omp(1,icore))
            dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + dbr * dlvec3_omp(1:dim1,icore)
            call lcpams_omp (iorbs2, ispin, isgn, 1, 2, icore)
            call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                              lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec4_omp(1,icore))
            dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) - dbr * dlvec4_omp(1:dim1,icore)
            if (dabs(dbi) .ge. dnano) then
              dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + dbi * dlvec4_omp(1:dim1,icore)
              dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) + dbi * dlvec3_omp(1:dim1,icore) 
            end if
!
            call lcpams_omp (iorbs2, ispin, isgn, 2, 2, icore)
            call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                              lentmp1_omp(icore), dlvec1_omp(1,icore), dlvec3_omp(1,icore))
            dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + ddr * dlvec3_omp(1:dim1,icore)
            call lcpams_omp (iorbs2, ispin, isgn, 2, 2, icore)
            call ams_dcoogemv('t', dim1, lvtmp1_omp(1,icore), lrtmp1_omp(1,icore), lctmp1_omp(1,icore), &
                              lentmp1_omp(icore), dlvec2_omp(1,icore), dlvec4_omp(1,icore))
            dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) - ddr * dlvec4_omp(1:dim1,icore)
            if (dabs(ddi) .ge. dnano) then
              dlvec7_omp(1:dim1,icore) = dlvec7_omp(1:dim1,icore) + ddi * dlvec4_omp(1:dim1,icore)
              dlvec8_omp(1:dim1,icore) = dlvec8_omp(1:dim1,icore) + ddi * dlvec3_omp(1:dim1,icore) 
            end if
          end if
        end do
      end if
!
      call l2hvec(nrho, dmtmp5_omp(1,1,icore), nrho, dlvec7_omp(1,icore))
      call l2hvec(nrho, dmtmp6_omp(1,1,icore), nrho, dlvec8_omp(1,icore))
!
!here, different threads may access the same lnj at the same time, and 
!thus operation on rcgrhs should be done one at a time by using 'omp critical'
!$omp critical
      rcgrhs(1:nrho, 1:nrho, lnj) = rcgrhs(1:nrho, 1:nrho, lnj) +                         &
                     dcmplx(dmtmp5_omp(1:nrho, 1:nrho, icore), dmtmp6_omp(1:nrho, 1:nrho, icore))
!$omp end critical
!      rcgrhs_omp(1:nrho,1:nrho,lnj,icore) = rcgrhs_omp(1:nrho,1:nrho,lnj,icore) +          & 
!                     dcmplx(dmtmp5_omp(1:nrho, 1:nrho, icore), dmtmp6_omp(1:nrho, 1:nrho, icore))
    end do
!
! nearest higher tier contribution (N + 1)
!
    if (itier .lt. ntier) then
!                 
      do isgn=1,nsgn
        jsgn = nsgn + 1 - isgn
        do ispin=1,nspin
          do iorbs=1,norbs
            do ialf=1,nalf
               do imats=1,ncor
                 lnj        = indexcoef(lnl)
                 itype      = itypecoef(lnl)
                 ifactfront = ifactfcoef(lnl)
                 ifactrear  = ifactrcoef(lnl)
                 lnl        = lnl + 1
!cycle ! +
                 dtmp1      = dble(ifactfront)
                 dtmp2      = dble(ifactrear )
                 if (lscale) then
                    dtmp1 = dtmp1 * dbsqrt(iorbs,ispin,imats,ialf,isgn)
                    dtmp2 = dtmp2 * dbsqrt(iorbs,ispin,imats,ialf,isgn)
                 end if
                 call liouamsmat_omp(1, itype, iorbs, ispin, jsgn,  dtmp1, rcgtmp(1,1,lni), nrho, &
                                     czero,  cmtmp2_omp(1,1,icore), nrho, icore)
                 call liouamsmat_omp(2, itype, iorbs, ispin, jsgn, -dtmp2, rcgtmp(1,1,lni), nrho, &
                                     cunity, cmtmp2_omp(1,1,icore), nrho, icore)
!$omp critical
                 rcgrhs(1:nrho, 1:nrho, lnj) = rcgrhs(1:nrho, 1:nrho, lnj) + cmtmp2_omp(1:nrho, 1:nrho, icore)
!$omp end critical
!                 rcgrhs_omp(1:nrho,1:nrho,lnj,icore) = rcgrhs_omp(1:nrho,1:nrho,lnj,icore) + &
!                                                       cmtmp2_omp(1:nrho, 1:nrho, icore)
               end do
            end do
          end do
        end do
      end do
    end if
    60 continue
  end do
!$omp end do nowait
!$omp end parallel
!
  call cpu_time(cpu2)
!  write(6,100)itier, nlast(itier), cpu2 - cpu1
!  call flush(6)
100 format('itier = ', I3, ' length ', I8, 2x, f12.6)
end do 
!
!do icore=1,nthreads_omp
!   rcgrhs(1:nrho,1:nrho,1:nunk) = rcgrhs(1:nrho,1:nrho,1:nunk) + rcgrhs_omp(1:nrho,1:nrho,1:nunk,icore)
!end do
!
! replace the first equation by normalization condition
!
! (1) remove the contribution to the equations corresponding to the
!     time derivative of rho(1,1,1)
!
dmtmp1(1:nrho, 1:nrho) =  dble(rcgtmp(1:nrho, 1:nrho, 1))
dmtmp2(1:nrho, 1:nrho) = dimag(rcgtmp(1:nrho, 1:nrho, 1))
call h2lvec(nrho, dmtmp1, nrho, dlvec1)
call h2lvec(nrho, dmtmp2, nrho, dlvec2)
!
dlmat1 = lhsrnl - lhsrnr
call modifylmat(nrho, dlmat1, dim1, 2)
dlmat1(2:dim1, 1:dim1) = 0.d0
dlmat3 = dlmat1
call modifylmat(nrho, dlmat1, dim1, 0)
call modifylmat(nrho, dlmat3, dim1, 1)
call dgemv('t', dim1, dim1, 1.d0, dlmat1, dim1, dlvec1, 1, 0.d0, dlvec3, 1)
call dgemv('t', dim1, dim1, 1.d0, dlmat3, dim1, dlvec2, 1, 0.d0, dlvec4, 1)
if (lhsi) then
  dlmat2 = lhsinl - lhsinr
  call modifylmat(nrho, dlmat2, dim1, 2)
  dlmat2(2:dim1, 1:dim1) = 0.d0
  dlmat4 = dlmat2
  call modifylmat(nrho, dlmat2, dim1, 0)
  call modifylmat(nrho, dlmat4, dim1, 1)
  call dgemv('t', dim1, dim1, 1.d0, dlmat2, dim1, dlvec2, 1, 1.d0, dlvec3, 1)
  call dgemv('t', dim1, dim1,-1.d0, dlmat4, dim1, dlvec1, 1, 1.d0, dlvec4, 1)
end if
call l2hvec(nrho, dmtmp5, nrho, dlvec3)
call l2hvec(nrho, dmtmp6, nrho, dlvec4)
rcgrhs(1:nrho, 1:nrho, 1) = rcgrhs(1:nrho, 1:nrho, 1) -                            &
                            dcmplx(dmtmp5(1:nrho, 1:nrho), dmtmp6(1:nrho, 1:nrho))
!
do lnj=ifirst(2), ilast(2)
  nj = indextable(lnj)
  if (mpm(nj) .ne. 1) then
    write(6,*)
    write(6,*)' error! wrong mpm for itier=2 ', lni, nj, mpm(nj)
    stop
  end if
!----
!cycle
!----  
  isgn   = 1
  jsgn   = 2
  ialf   = ilead(nj)
  iorbs  = morbs(nj)
  ispin  = mspin(nj)
  imats  = mcor(nj)
!
  if (lscale) then
     dtmp1 = dbsqrt(iorbs,ispin,imats,ialf,isgn)
     dtmp2 = dbsqrt(iorbs,ispin,imats,ialf,jsgn)
  else
     dtmp1 = 1.d0
     dtmp2 = 1.d0
  end if
!
  call look4drawer(isgn, ialf, iorbs, ispin, imats, nn)
  lni = nfirst(2) - 1 + nn
  dlvec7 = 0.d0   ! Yr
  dlvec8 = 0.d0   ! Yi
!
  call lcpams (iorbs, ispin, 2, 1, 1)
  call modlams(iorbs, ispin, 2, 1, 1, 2)
  do nk=1,lentmp1
    if (lrtmp1(nk) .ge. 2) lvtmp1(nk) = 0.d0
  end do
  call ams_dcoogemv('t', dim1, lvtmp1(1), lrtmp1(1), lctmp1(1), lentmp1, dlvec2, dlvec3)
  call ams_dcoogemv('t', dim1, lvtmp1(1), lrtmp1(1), lctmp1(1), lentmp1, dlvec1, dlvec4)
  dlvec7 = dlvec7 - dlvec3 * dtmp1
  dlvec8 = dlvec8 + dlvec4 * dtmp1
!
  call lcpams (iorbs, ispin, 2, 2, 1)
  call modlams(iorbs, ispin, 2, 2, 1, 2)
  do nk=1,lentmp1
    if (lrtmp1(nk) .ge. 2) lvtmp1(nk) = 0.d0
  end do
  call ams_dcoogemv('t', dim1, lvtmp1(1), lrtmp1(1), lctmp1(1), lentmp1, dlvec2, dlvec3)
  call ams_dcoogemv('t', dim1, lvtmp1(1), lrtmp1(1), lctmp1(1), lentmp1, dlvec1, dlvec4)
  dlvec7 = dlvec7 + dlvec3 * dtmp1
  dlvec8 = dlvec8 - dlvec4 * dtmp1
!
  call lcpams (iorbs, ispin, 1, 1, 2)
  call modlams(iorbs, ispin, 1, 1, 2, 2)
  do nk=1,lentmp1
    if (lrtmp1(nk) .ge. 2) lvtmp1(nk) = 0.d0
  end do
  call ams_dcoogemv('t', dim1, lvtmp1(1), lrtmp1(1), lctmp1(1), lentmp1, dlvec2, dlvec3)
  call ams_dcoogemv('t', dim1, lvtmp1(1), lrtmp1(1), lctmp1(1), lentmp1, dlvec1, dlvec4)
  dlvec7 = dlvec7 - dlvec3 * dtmp2
  dlvec8 = dlvec8 - dlvec4 * dtmp2
!
  call lcpams (iorbs, ispin, 1, 2, 2)
  call modlams(iorbs, ispin, 1, 2, 2, 2)
  do nk=1,lentmp1
    if (lrtmp1(nk) .ge. 2) lvtmp1(nk) = 0.d0
  end do
  call ams_dcoogemv('t', dim1, lvtmp1(1), lrtmp1(1), lctmp1(1), lentmp1, dlvec2, dlvec3)
  call ams_dcoogemv('t', dim1, lvtmp1(1), lrtmp1(1), lctmp1(1), lentmp1, dlvec1, dlvec4)
  dlvec7 = dlvec7 + dlvec3 * dtmp2
  dlvec8 = dlvec8 + dlvec4 * dtmp2
!
  call l2hvec(nrho, dmtmp5, nrho, dlvec7)
  call l2hvec(nrho, dmtmp6, nrho, dlvec8)
  rcgrhs(1:nrho, 1:nrho, lni) = rcgrhs(1:nrho, 1:nrho, lni) -                         &
                                dcmplx(dmtmp5(1:nrho, 1:nrho), dmtmp6(1:nrho, 1:nrho)) 
end do
!
22 continue
!
! (2) impose normalization condition
!
do ni=1,nrho
  rcgrhs(ni,ni,1) = rcgrhs(ni,ni,1) + rcgtmp(1,1,1)
end do
!
300 continue
!
do ni=1,nrho
  do nj=ni+1,nrho
    rcgrhs(ni,nj,1) = czero
  end do
end do
do ni=1,nrho
  rcgrhs(ni,ni,1) = dcmplx(dble(rcgrhs(ni,ni,1)), 0.d0)
end do
!
end subroutine calc_atrx_bicg_omp

