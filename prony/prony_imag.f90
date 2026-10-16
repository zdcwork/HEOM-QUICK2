subroutine prony_job_imag(sigma_tmp, prony_tmp, dscaling)
  !
  implicit none
  include "../include/sizes"
  include "../include/common"
  !
  integer, intent(inout) :: prony_tmp
  real*8, intent(inout)  :: sigma_tmp, dscaling
  !
  integer    :: ni, nj, nk, nl, nm, ncount, nsign
  integer    :: dimension_Hankel, tt_sample_for_prony
  integer    :: ntmp1, ntmp2, ntmp3
  integer, allocatable    :: nvectmp1(:)
  integer, allocatable    :: nindextmp(:)
  !
  real*8     :: sigma
  real*8     :: tt_min_for_prony, tt_step_for_prony
  real*8     :: tt_max_for_prony
  real*8     :: dtmp1, dtmp2, dtmp3, error
  real*8, allocatable    :: tt_for_prony(:)
  real*8, allocatable :: dvectmp1(:), dvectmp2(:)
  !
  integer  :: nprony, istat, info, lwork, liwork
  integer, allocatable  :: iwork(:)
  complex*16 :: ctmp1, ctmp2, ctmp3, cunit
  real*8, allocatable    :: hankel(:,:)
  real*8, allocatable    :: hankel_tmp_eigval(:), hankel_tmp_eigvec(:,:)
  real*8, allocatable    :: dwork(:)
  complex*16, allocatable :: cpara(:), cpara_mat(:,:)
  complex*16, allocatable :: eta_prony(:), diss_rate_prony(:), res_prony(:)
  complex*16, allocatable :: cvectmp1(:), cvectmp2(:), cdwork(:)
  complex*16, allocatable :: cmattmp1(:,:), cmattmp2(:,:), cmattmp3(:,:), cmattmp4(:,:), cmattmp5(:,:)
  complex*16, allocatable :: zhankel_eigval(:), zhankel_eigvec(:,:)
  complex*16, allocatable :: prony_roots(:)
  logical :: lfile_exist
  character(len=128) :: buffer
  
  !!!!!!!!!!!! necessary and important constant !!!!!!!!!!!!
  
  nsign = 2
  tt_min_for_prony = 0.0d0
  cunit = dcmplx(1.d0, 0.d0)
  
  !!!!!!!!!!!! input parameters !!!!!!!!!!!!
 
  sigma = sigma_tmp
  nprony = prony_tmp
  !dimension_Hankel = 3250       !!! default value = 1500
  !tt_max_for_prony = 350.0d0    !!! tt_max_for_prony should be sufficiently large to describe long time dynamics
  
  !!!!!!!!!!!! build Hankel matirx for C(t) !!!!!!!!!!!!
 
  write(6,*) 'prony_job: reading data from input file' 
  inquire(file='IFT.data',exist=lfile_exist)
  if (lfile_exist) then
    write(6,*) 'prony_job: reading data from IFT.data'
    open(unit=13, file='IFT.data')
    rewind(13)
    tt_sample_for_prony = -1
    do while(.not. is_iostat_end(istat))
      read(13,'(a)', iostat=istat) buffer
      tt_sample_for_prony = tt_sample_for_prony + 1
    end do
    write(6,*) 'prony_job: the size of input data', tt_sample_for_prony
    dimension_Hankel = int(tt_sample_for_prony + 1) / 2
    write(6,*) 'prony_job: dimension of Hankel matrix', dimension_Hankel
    rewind(13)
        !
    allocate(tt_for_prony(tt_sample_for_prony),&
             res_prony(tt_sample_for_prony), STAT=istat)
    res_prony = czero
    tt_for_prony = 0.d0
    do ni=1, tt_sample_for_prony
      read(13,200) dtmp1,dtmp2,dtmp3
      tt_for_prony(ni) = dtmp1
      res_prony(ni) = dcmplx(dtmp2,dtmp3) * dscaling
    end do
    !
  end if
  200 format(3(E23.16,1x))
  !
  write(6,*) 'prony_job: initial time =', tt_for_prony(1)
  write(6,*) 'prony_job: time step   =', tt_for_prony(2)
  write(6,*) 'prony_job: final time  =', tt_for_prony(tt_sample_for_prony)
  !write(6,*)
  
  allocate(hankel(dimension_Hankel, dimension_Hankel), STAT=istat)
  hankel = 0.d0
  do ni = 1, dimension_Hankel
    do nj = 1, dimension_Hankel
      hankel(ni,nj) = dimag(res_prony(ni + nj - 1))
    end do
  end do
   
  ! add the symmetry check at here!
  
  !!!!!!!!!!!! Takagi factorization for input Hankel matirx !!!!!!!!!!!!
  
  allocate(hankel_tmp_eigval(dimension_Hankel), &
           hankel_tmp_eigvec(dimension_Hankel,dimension_Hankel), STAT=istat)
  hankel_tmp_eigvec = hankel
  hankel_tmp_eigval = 0.d0
  
  lwork = (1 + 6 * dimension_Hankel + 2 * dimension_Hankel**2) * 2
  liwork = (3 + 5 * dimension_Hankel) * 2
  allocate(dwork(lwork), iwork(liwork), STAT=istat)
  dwork = 0.0d0
  iwork = 0
  !call dsyevd(JOBZ, UPLO, N, A( lda, N ), LDA, W( N ), WORK( LWORK ), LWORK, IWORK( LIWORK ), LIWORK, NFO)
  call dsyevd('V','U', dimension_Hankel, hankel_tmp_eigvec,&
              dimension_Hankel, hankel_tmp_eigval,&
              dwork, lwork, iwork, liwork, info)
  if (info .ne. 0) then
   write(6,*)'prony_job: error! eigenvalue calculation failed for hankel matrix', info
   stop
  end if
  
  allocate(nindextmp(dimension_Hankel), STAT=istat)
  allocate(dvectmp1(dimension_Hankel), STAT=istat)
  allocate(cmattmp1(dimension_Hankel,dimension_Hankel),&
           cmattmp2(dimension_Hankel,dimension_Hankel), STAT=istat)
  allocate(zhankel_eigval(dimension_Hankel), &
           zhankel_eigvec(dimension_Hankel,dimension_Hankel), STAT=istat)
  
  nindextmp = 0
  dvectmp1 = 0.0d0
  cmattmp1 = czero
  cmattmp2 = czero
  
  do nj = 1, dimension_Hankel
    if (hankel_tmp_eigval(nj) < 0.0d0) then
      do ni = 1, dimension_Hankel
        zhankel_eigvec(ni,nj) = -eye * hankel_tmp_eigvec(ni,nj)
      end do
    else if (hankel_tmp_eigval(nj) > 0.0d0) then
      do ni = 1, dimension_Hankel
        zhankel_eigvec(ni,nj) = cunit * hankel_tmp_eigvec(ni,nj)
      end do
    end if
  end do
  
  do ni = 1, dimension_Hankel
    dvectmp1(ni) = dabs(hankel_tmp_eigval(ni))
    nindextmp(ni) = ni
  end do
 
  
  deallocate(hankel_tmp_eigval, hankel_tmp_eigvec, STAT=istat)
  
  cmattmp1 = zhankel_eigvec
  
  ! sort dabs(hankel_eigval) in the descending order
  do ni = 1, dimension_Hankel -1
    do nj = ni, 1, -1
      if (dvectmp1(nj) < dvectmp1(nj + 1)) then
        dtmp1 = dvectmp1(nj + 1)
        dvectmp1(nj + 1) = dvectmp1(nj)
        dvectmp1(nj) = dtmp1
        !
        ntmp1 = nindextmp(nj + 1)
        nindextmp(nj + 1) = nindextmp(nj)
        nindextmp(nj) = ntmp1
      end if
     end do
  end do
  
  do ni = 1, dimension_Hankel
    do nj = 1, dimension_Hankel
      cmattmp2(ni,nj) = cmattmp1(ni,nindextmp(nj))
    end do
  end do
  
  zhankel_eigval = dcmplx(dvectmp1, 0.0d0)
  zhankel_eigvec = cmattmp2
 
  
  deallocate(nindextmp, STAT=istat)
  deallocate(dvectmp1, STAT=istat)
  deallocate(cmattmp1, cmattmp2, STAT=istat)
  
  !!!!!!!!!!!! check the accurancy of Takagi factorization for imaginary Hankel matirx !!!!!!!!!!!!
 
  allocate(cmattmp1(dimension_Hankel,dimension_Hankel),&
           cmattmp2(dimension_Hankel,dimension_Hankel), STAT=istat)
  allocate(cmattmp3(dimension_Hankel,dimension_Hankel),&
           cmattmp4(dimension_Hankel,dimension_Hankel), STAT=istat)
  cmattmp1 = czero
  cmattmp2 = czero
  cmattmp3 = czero
  cmattmp4 = czero
  error = 0.0d0
  
  ! cmattmp1 = diagional matirx (hankel_eigval)
  do ni = 1, dimension_Hankel
    cmattmp1(ni,ni) = zhankel_eigval(ni)
    do nj = 1, dimension_Hankel
      cmattmp4(ni,nj) = zhankel_eigvec(nj,ni)
    end do
  end do
  
  ! zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  ! C = alpha * A * B + beta * C
  ! cmattmp2 =  zhankel_eigvec * hankel_eigval 
  call zgemm('N', 'T', dimension_Hankel, dimension_Hankel,&
             dimension_Hankel, cunit, zhankel_eigvec, &
             dimension_Hankel, cmattmp1, dimension_Hankel,&
             czero, cmattmp2, dimension_Hankel)
        
   
  ! cmattmp3 =  zhankel_eigvec * hankel_eigval * zhankel_eigvec.T
  call zgemm('N', 'T', dimension_Hankel, dimension_Hankel,&
             dimension_Hankel, cunit, cmattmp2, &
             dimension_Hankel, zhankel_eigvec, dimension_Hankel,&
             czero, cmattmp3, dimension_Hankel)
  
  dtmp1 = 0.d0
  do ni = 1, dimension_Hankel
    do nj = 1, dimension_Hankel
      dtmp1 = max(dtmp1, cdabs( dcmplx(hankel(ni,nj), 0.0d0) - cmattmp3(ni,nj) ))
      error = error + cdabs( dcmplx(hankel(ni,nj), 0.0d0) - cmattmp3(ni,nj) )
    end do
  end do
  
  !write(6,*)'prony_job: error of Takagi facotization for imagnary part', error 
  !write(6,*)'prony_job: the maxium error element is', dtmp1
  !write(6,*)
  
  deallocate(cmattmp1,cmattmp2,cmattmp3,cmattmp4, STAT=istat)
  
  !!!!!!!!!!!! calculate the roots of polynomial !!!!!!!!!!!!
  
  allocate(cpara(dimension_Hankel - 1), &
           cpara_mat(dimension_Hankel - 1,dimension_Hankel - 1), STAT=istat)
  allocate(prony_roots(dimension_Hankel - 1), STAT=istat)
  cpara = czero
  cpara_mat = czero
  prony_roots = czero
  
  do ni = 1, dimension_Hankel - 1
    cpara(ni) = zhankel_eigvec(ni, nprony + 1) / &
                (zhankel_eigvec(dimension_Hankel, nprony + 1))
    !if (ni <= 10 .or. (dimension_Hankel - 1 - ni) <= 10) then
    !  write(6,*)'numerator', zhankel_eigvec(ni, nprony + 1)
    !  write(6,*)'denominator',(zhankel_eigvec(dimension_Hankel, nprony + 1))
    !end if
  end do
  
  
  do ni = 1, dimension_Hankel - 1
    do nj = 1, dimension_Hankel - 1
      if (nj == (ni + 1)) then
        cpara_mat(ni,nj) = cunit
      end if
    end do
    if (ni == dimension_Hankel - 1) then
      do nj = 1, dimension_Hankel - 1
        cpara_mat(ni,nj) = -cunit * cpara(nj)
      end do
    end if
  end do
  !
  lwork = 5*(dimension_Hankel-1)
  allocate(cmattmp1(dimension_Hankel - 1, dimension_Hankel - 1),&
           cmattmp2(dimension_Hankel - 1, dimension_Hankel - 1), STAT=istat)
  allocate(cdwork(lwork), dwork(2*(dimension_Hankel-1)), STAT=istat)
  cdwork = czero
  cmattmp1 = czero
  cmattmp2 = czero
  dwork = 0.0d0
  dtmp1 = 0.0d0
  !
  !call zgeev(JOBVL, JOBVR, N, A( lda, N ), LDA, W( N ), VL( ldvl, N ), LDVL, VR( ldvr, N ), LDVR, WORK( LWORK ), LWORK, RWORK( 2N ), INFO)
  call zgeev('N','N', dimension_Hankel - 1, cpara_mat,&
             dimension_Hankel - 1, prony_roots, &
             cmattmp1, dimension_Hankel - 1, cmattmp2,& 
             dimension_Hankel - 1, cdwork, 5*(dimension_Hankel-1), dwork, info)
  !
  deallocate(cdwork, STAT=istat)
  deallocate(cmattmp1,cmattmp2, STAT=istat)
  deallocate(zhankel_eigvec, zhankel_eigval, STAT=istat)
  
  !!!!!!!!!!!! get the number of nonzero roots for imagniary part !!!!!!!!!!!!
  
  allocate(cvectmp1(dimension_Hankel), STAT=istat)
  cvectmp1 = czero
  ncount = 0
  do ni = 1, dimension_Hankel - 1
    if (cdabs(prony_roots(ni)) >= 10E-16 ) then
      ncount = ncount + 1
      cvectmp1(ncount) = prony_roots(ni)
    end if
  end do
  !
  if (ncount == (dimension_Hankel - 1)) then
   !write(6,*) 'prony_job: there is no zero element'
   !write(6,*) 'prony_job: number of nonzero element', ncount
   !write(6,*)
  else if (ncount .ne. (dimension_Hankel - 1)) then
    deallocate(prony_roots)
    allocate(prony_roots(ncount), STAT=istat)
    !write(6,*) 'prony_job: number of nonzero element', ncount
    !write(6,*)
    do ni = 1, ncount
      prony_roots(ni) = cvectmp1(ni)
    end do
  end if
  
  deallocate(cvectmp1, STAT=istat)
      
  !!!!!!!!!!!! get the order of dabs(prony_roots) in dabs(prony_roots) sorted in the ascending order !!!!!!!!!!!! 
     
  allocate(nindextmp(ncount), STAT=istat)
  allocate(dvectmp1(ncount), dvectmp2(ncount), STAT=istat)
  dvectmp1 = 0.0d0
  dvectmp2 = 0.0d0
  nindextmp = 0
  do ni = 1, ncount
    dvectmp1(ni) = cdabs(prony_roots(ni))
    nindextmp(ni) = ni
    !re_nindextmp(ni) = ni
  end do
  
  dvectmp2 = dvectmp1
  do ni = ncount -1, 1, -1
    do nj = 1, ni
      if (dvectmp2(nj) > dvectmp2(nj + 1)) then
        dtmp1 = dvectmp2(nj)
        dvectmp2(nj) = dvectmp2(nj + 1)
        dvectmp2(nj + 1) = dtmp1
        !
        ntmp1 = nindextmp(nj)
        nindextmp(nj) = nindextmp(nj + 1)
        nindextmp(nj + 1) = ntmp1
      end if
    end do
  end do
  
  allocate(cvectmp1(nprony), STAT=istat)
  cvectmp1 = czero
  
  do ni = 1, nprony
    cvectmp1(ni) = prony_roots(nindextmp(ni))
  end do
  
  deallocate(prony_roots, STAT=istat)
  
  !!!!!!!!!!!! calculate the dissipation rates of imagniary part !!!!!!!!!!!!
  
  allocate(prony_roots(nprony), diss_rate_prony(nprony), STAT=istat)
  prony_roots = cvectmp1
  diss_rate_prony = czero
  !
  cvectmp1 = czero
  do ni = 1, nprony
    ctmp1 = dcmplx(atan( dimag(prony_roots(ni)), dble(prony_roots(ni))), 0.0d0)
    cvectmp1(ni) = -2.d0 * (dimension_Hankel - 1) / tt_for_prony(tt_sample_for_prony) &
              * ( dlog( cdabs(prony_roots(ni)) ) + eye * ctmp1)
  end do
  
  diss_rate_prony = cvectmp1
  deallocate(cvectmp1, STAT=istat)
  
  !!!!!!!!!!!! calculate the coefficients of imagniary part !!!!!!!!!!!!
  ! calculate the coefficients by the least squares method
  
  allocate(cmattmp1(tt_sample_for_prony, nprony), STAT=istat)
  allocate(cmattmp2(nprony, tt_sample_for_prony), STAT=istat)
  allocate(cmattmp3(nprony, nprony), STAT=istat)
  allocate(cmattmp4(nprony, nprony), STAT=istat)
  allocate(cvectmp1(tt_sample_for_prony), STAT=istat)
  allocate(nvectmp1(nprony), STAT=istat)
  cmattmp1 = czero
  cmattmp2 = czero
  cmattmp3 = czero
  cmattmp4 = czero
  cvectmp1 = czero
  nvectmp1 = 0
  
  do ni = 1, tt_sample_for_prony
    do nj = 1, nprony
      cmattmp1(ni,nj) = prony_roots(nj)**(ni - 1)
      cmattmp2(nj,ni) = prony_roots(nj)**(ni - 1)
    end do
  end do
  
  ! cmattmp3 = cmattmp2 * cmattmp1 = omega.T * omega
  call zgemm('N', 'N', nprony, nprony, tt_sample_for_prony, cunit, cmattmp2, &
           nprony, cmattmp1, tt_sample_for_prony, czero, cmattmp3, nprony)
  
  ! cvectmp1 = cmattmp2 * dimag(res_prony) = omega.T * dimag(res_prony)
  do ni = 1, nprony
    do nj = 1, tt_sample_for_prony
      cvectmp1(ni) = cvectmp1(ni) + cmattmp2(ni,nj) * dimag(res_prony(nj))
    end do
  end do
  
  ! inverse cmattmp3
  lwork = 5 * nprony
  allocate(cvectmp2(lwork), STAT=istat)
  cvectmp2 = czero
  cmattmp4 = cmattmp3
  
  ! calculate the coefficients by the least squares method
  call zgetrf(nprony, nprony, cmattmp3, nprony, nvectmp1, info)
  call zgetri(nprony, cmattmp3, nprony, nvectmp1, cvectmp2, lwork, info) 
  
  deallocate(nvectmp1, cvectmp2, STAT=istat)
  
  
  allocate(cmattmp5(nprony, nprony), STAT=istat)
  
  cmattmp5 = czero
  call zgemm('N', 'N', nprony, nprony, nprony, cunit, cmattmp4, &
           nprony, cmattmp3, nprony, czero, cmattmp5, nprony)
  
  allocate(eta_prony(nprony), STAT=istat)
  eta_prony = czero
  
  ! cmattmp3 * cvectmp1
  do ni = 1, nprony
    do nj = 1, nprony
      eta_prony(ni) = eta_prony(ni) + cmattmp3(ni,nj) * cvectmp1(nj)
    end do
  end do
  eta_prony = eta_prony * eye
  eta_prony = eta_prony / dscaling
  
  deallocate(cmattmp1, cmattmp2, cmattmp3,&
             cmattmp4, cmattmp5, cvectmp1, prony_roots, STAT=istat)
  
  !!!!!!!!!!!! output dissipation rate and coefficients !!!!!!!!!!!!
 
  !do ni = 1, nprony
  !  dtmp1 = dble(eta_prony(ni))
  !  dtmp2 = dimag(eta_prony(ni))
  !  if (dabs(dtmp1) <= 1.d-6) then
  !      dtmp1 = 0.d0
  !  end if
  !  if (dabs(dtmp2) <= 1.d-6) then
  !      dtmp2 = 0.d0
  !  end if
  !  eta_prony(ni) = dcmplx(dtmp1,dtmp2)
  !  
  !  dtmp1 = dble(diss_rate_prony(ni))
  !  dtmp2 = dimag(diss_rate_prony(ni))
  !  if (dabs(dtmp1) <= 1.d-6) then
  !      dtmp1 = 0.d0
  !  end if
  !  if (dabs(dtmp2) <= 1.d-6) then
  !      dtmp2 = 0.d0
  !  end if
  !  diss_rate_prony(ni) = dcmplx(dtmp1,dtmp2)
  !end do
  
 open(unit = 94, file='PRONY_CB.data', POSITION='APPEND')
 !rewind(94)
 !write(94,*)'icor, Re[cb], Imag[cb]'
 open(unit = 95, file='PRONY_CD.data', POSITION='APPEND')
 !rewind(95)
 !write(95,*)'icor, Re[cd], Imag[cd]'
 open(unit = 96, file='PRONY_GAMA.data', POSITION='APPEND')
 !rewind(96)
 !write(96,*)'icor, Re[cgama], Imag[cgama]'
 do nk = 1, nprony
   write(94,1001)dble(eta_prony(nk))/hbar**2,&
    dimag(eta_prony(nk))/hbar**2
   write(95,1001)dble(eta_prony(nk))/hbar**2,&
    -dimag(eta_prony(nk))/hbar**2
   write(96,1001)-dble(diss_rate_prony(nk))/hbar,&
    -dimag(diss_rate_prony(nk))/hbar
 end do
 close(94)
 close(95)
 close(96)
 
 1001 format(2(E23.16,1x))
end subroutine prony_job_imag
