subroutine prony_method(sigma_tmp,prony_real_tmp,prony_imag_tmp)
 !
 implicit none
 !
 include "../include/sizes"
 include "../include/common"
 !
 integer, intent(inout) :: prony_real_tmp, prony_imag_tmp
 real*8, intent(inout)  :: sigma_tmp
 integer :: prony_real, prony_imag
 !
 integer    :: ni, nj, nk, nl, nm, ncount, nsign
 integer    :: npade 
 integer    :: dimension_Hankel, tt_sample_for_prony
 integer    :: ntmp1, ntmp2, ntmp3
 integer, allocatable    :: nvectmp1(:)
 integer, allocatable    :: nindextmp(:)
 
 !real*8     :: hbar, pi
 real*8     :: sys_bath_coupling, band_width, band_center
 real*8     :: chemical_potential, T, sigma
 real*8     :: tt_min_for_prony, tt_step_for_prony
 real*8     :: tt_max_for_prony
 real*8     :: dtmp1, dtmp2, dtmp3, error
 real*8, allocatable    :: tt_for_prony(:)
 real*8, allocatable :: dvectmp1(:), dvectmp2(:)
 
 logical  :: lpade, lfile_exist, lonly_real, lonly_imag
 character(len=128) :: buffer
 integer  :: istat, info, lwork, liwork
 integer, allocatable  :: iwork(:)
 !complex*16 :: cunit, eye, czero
 complex*16 :: ctmp1, ctmp2, ctmp3
 real*8, allocatable    :: hankel_real(:,:), hankel_imag(:,:)
 real*8, allocatable    :: hankel_tmp_eigval(:), hankel_tmp_eigvec(:,:)
 real*8, allocatable    :: dwork(:)
 complex*16, allocatable :: zpeigv(:), zpcoef(:), zppole_tmp(:,:)
 complex*16, allocatable :: zppole(:), cpara(:), cpara_mat(:,:)
 complex*16, allocatable :: eta_psd(:), eta_psd_dag(:), diss_rate_psd(:)
 complex*16, allocatable :: eta_prony(:), diss_rate_prony(:)
 complex*16, allocatable :: res_prony(:)
 complex*16, allocatable :: cvectmp1(:), cvectmp2(:), cdwork(:)
 complex*16, allocatable :: cmattmp1(:,:), cmattmp2(:,:), cmattmp3(:,:)
 complex*16, allocatable :: cmattmp4(:,:), cmattmp5(:,:)
 complex*16, allocatable :: zhankel_real_eigval(:)
 complex*16, allocatable :: zhankel_imag_eigval(:)
 complex*16, allocatable :: zhankel_real_eigvec(:,:)
 complex*16, allocatable :: zhankel_imag_eigvec(:,:)
 complex*16, allocatable :: prony_roots_real(:), prony_roots_imag(:)
 complex*16, allocatable :: diss_rate_prony_real(:)
 complex*16, allocatable :: diss_rate_prony_imag(:)
 complex*16, allocatable :: eta_prony_real(:), eta_prony_imag(:)
 complex*16, allocatable :: eta_prony_tmp(:), diss_rate_prony_tmp(:)
 
 !!!!!!!!!!!! necessary and important constant !!!!!!!!!!!!
 
 !hbar = 0.658211928d0
 !pi = 3.14159265358979323846d0
 !cunit = (1.0d0, 0.0d0)
 !eye = (0.0d0, 1.0d0)
 !czero = (0.0d0, 0.0d0)
 lpade = .False.
 lonly_real = .False.
 lonly_imag = .False.
 
 nsign = 2
 tt_min_for_prony = 0.0d0
 
 !!!!!!!!!!!! input parameters !!!!!!!!!!!!
 
 npade = 650
 sigma = int(sigma_tmp)
 
 prony_real = prony_real_tmp
 prony_imag = prony_imag_tmp
 !dimension_Hankel = 3225      ! default value = 1500
 !tt_max_for_prony = 184.0d0   ! default value = 80.0d0 
                              ! Notice that tt_max_for_prony should
                              ! not be too large
 !sys_bath_coupling = 1.0d0    ! eV, linewidth(norbs,nalf) in HEOM-QUICK
 !band_width = 5.0d0           ! eV, gammap(nvar2,nspin,malf) 
                              ! or bandwidth(nalf) in HEOM-QUICK
 !band_center = 0.d0           ! eV, dwp(iorbs2,ispin,ialf) 
                              ! or bcenter(nalf,nspin) in HEOM-QUICK
 !T = 0.001d0                  ! eV, dinvbeta(nalf) in HEOM-QUICK
 
 !open(unit=5,file='input',status='unknown')
 !rewind(5)
 !read(UNIT=5,NML=prony,END=50)
 !50 continue
 !close(5)
 !
 !!!!!!!!!!!! get poles and coefficients !!!!!!!!!!!!
 if(lpade) then
   write(6,*) 'prony_method: not compatible with Pade mode'
   goto 2000 
 else
   write(6,*) 'prony_method: reading data from input file' 
   inquire(file='IFT.data',exist=lfile_exist)
   if (lfile_exist) then
     write(6,*) 'prony_method: reading data from IFT.data'
     open(unit=13, file='IFT.data')
     rewind(13)
     read(13,'(a)') buffer
     tt_sample_for_prony = -1
     do while(.not. is_iostat_end(istat))
       read(13,'(a)', iostat=istat) buffer
       tt_sample_for_prony = tt_sample_for_prony + 1
     end do
     write(6,*) 'prony_method: the size of input data', tt_sample_for_prony
     if (mod(tt_sample_for_prony,2) .eq. 0) then
       write(6,*) 'prony_method: fatal error:',&
                 ' tt_sample_for_prony should be even'
       stop
     end if
     dimension_Hankel = (tt_sample_for_prony + 1) / 2
     write(6,*) 'prony_method: dimension of Hankel matrix', dimension_Hankel
     rewind(13)
         !
     allocate(tt_for_prony(tt_sample_for_prony),&
              res_prony(tt_sample_for_prony), STAT=istat)
     res_prony = czero
     tt_for_prony = 0.d0
     read(13,'(a)') buffer
     do ni=1, tt_sample_for_prony
       read(13,200) tt_for_prony(ni),dtmp1,dtmp2,dtmp3
       res_prony(ni) = dcmplx(dtmp1,dtmp2)
     end do
     !
     write(6,*) 'prony_method: prony_job: initial time =',&
                 tt_for_prony(1)
     write(6,*) 'prony_method: prony_job: time step   =',&
                 tt_for_prony(2) - tt_for_prony(1)
     write(6,*) 'prony_method: prony_job: final time  =',&
                 tt_for_prony(tt_sample_for_prony)
     write(6,*)
     !
   end if
 end if
 200 format(4(2x,E20.11))
 
 !!!!!!!!!!!! build real and imaginary Hankel matirx for C(t)
 !!!!!!!!!!!!!
 write(6,*) 'prony_method: building real and imaginary Hankel matrices'
 allocate(hankel_real(dimension_Hankel, dimension_Hankel),&
          hankel_imag(dimension_Hankel, dimension_Hankel), STAT=istat)
 hankel_real = 0.d0
 hankel_imag = 0.d0
 do ni = 1, dimension_Hankel
   do nj = 1, dimension_Hankel
     hankel_real(ni,nj) = dreal(res_prony(ni + nj - 1))
     hankel_imag(ni,nj) = dimag(res_prony(ni + nj - 1))
   end do
 end do
 write(6,*) 'prony_method: building real and imaginary Hankel matrices done'
 write(6,*)

 if (lonly_imag .and. .not. lonly_real) then
   prony_real = 0
   write(6,*) 'prony_method: only performing Prony for imag parts'
   !write(6,*) 'number of exp. func. for real part',&
   !            prony_real
   !write(6,*) 'number of exp. func. for imagniary part',&
   !            prony_imag
   !write(6,*)
 end if
 ! 
 if (lonly_real .and. .not. lonly_imag) then
   prony_imag = 0
   write(6,*) 'prony_method: only performing Prony for real parts'
   write(6,*) 'prony_method: number of exp. func. for real part',&
               prony_real
   write(6,*) 'prony_method: number of exp. func. for imagniary part',&
               prony_imag
   write(6,*)
   goto 1500
 end if
 !
 write(6,*) 'prony_method: number of exp. func. for real part',&
             prony_real
 write(6,*) 'prony_method: number of exp. func. for imagniary part',&
             prony_imag
 write(6,*)
 ! add the symmetry check at here!
 
 !!!!!!!!!!!! Takagi factorization for imaginary Hankel matirx
 !!!!!!!!!!!!!

 write(6,*) 'prony_method: performing Takagi factorization for imaginary part' 
 allocate(hankel_tmp_eigval(dimension_Hankel),&
          hankel_tmp_eigvec(dimension_Hankel,dimension_Hankel),&
          STAT=istat)
 hankel_tmp_eigvec = hankel_imag
 hankel_tmp_eigval = 0.d0
 
 lwork = (1 + 6 * dimension_Hankel + 2 * dimension_Hankel**2) * 2
 liwork = (3 + 5 * dimension_Hankel) * 2
 allocate(dwork(lwork), iwork(liwork), STAT=istat)
 dwork = 0.0d0
 iwork = 0
 !call dsyevd(JOBZ, UPLO, N, A( lda, N ), LDA, W( N ), WORK( LWORK ),
 !LWORK, IWORK( LIWORK ), LIWORK, NFO)
 write(6,*) 'prony_method: solving eigen-value of Hankel matrix'
 call dsyevd('V','U', dimension_Hankel, hankel_tmp_eigvec,&
             dimension_Hankel, hankel_tmp_eigval, dwork, lwork, iwork,&
             liwork, info)
 write(6,*) 'prony_method: eigen-value calculation done'
!!! DEBUG
 deallocate(dwork, iwork, STAT=istat)
!!! DEBUG
 if (info .ne. 0) then
   write(6,*)'prony_job: error!',& 
             ' eigenvalue calculation failed for hankel_imag',& 
               info
   stop
 end if
 !
 allocate(nindextmp(dimension_Hankel), STAT=istat)
 allocate(dvectmp1(dimension_Hankel), STAT=istat)
 allocate(cmattmp1(dimension_Hankel,dimension_Hankel),&
          cmattmp2(dimension_Hankel,dimension_Hankel), STAT=istat)
 allocate(zhankel_imag_eigval(dimension_Hankel),&
          zhankel_imag_eigvec(dimension_Hankel,dimension_Hankel),&
          STAT=istat)
 !
 nindextmp = 0
 dvectmp1 = 0.0d0
 cmattmp1 = czero
 cmattmp2 = czero
 
 do nj = 1, dimension_Hankel
   if (hankel_tmp_eigval(nj) < 0.0d0) then
     do ni = 1, dimension_Hankel
       zhankel_imag_eigvec(ni,nj) = -eye *&
       hankel_tmp_eigvec(ni,nj)
     end do
   else if (hankel_tmp_eigval(nj) > 0.0d0) then
     do ni = 1, dimension_Hankel
       zhankel_imag_eigvec(ni,nj) = cunity *&
       hankel_tmp_eigvec(ni,nj)
     end do
   end if
 end do
! 
 do ni = 1, dimension_Hankel
   dvectmp1(ni) = dabs(hankel_tmp_eigval(ni))
   nindextmp(ni) = ni
 end do
!
 deallocate(hankel_tmp_eigval, hankel_tmp_eigvec, STAT=istat)
 write(6,*) 'prony_method: Takagi factorization done for imaginary part' 
 !write(6,*)
! 
 cmattmp1 = zhankel_imag_eigvec
! 
! sort dabs(hankel_imag_eigval) in the descending order
 write(6,*) 'prony_method: sorting abs of eigen-values'
 do ni = 1, dimension_Hankel - 1
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
 !
 do ni = 1, dimension_Hankel
   do nj = 1, dimension_Hankel
     cmattmp2(ni,nj) = cmattmp1(ni,nindextmp(nj))
   end do
 end do
 !
 zhankel_imag_eigval = dcmplx(dvectmp1, 0.0d0)
 zhankel_imag_eigvec = cmattmp2
 !
 deallocate(nindextmp, STAT=istat)
 deallocate(dvectmp1, STAT=istat)
 deallocate(cmattmp1, cmattmp2, STAT=istat)
 write(6,*) 'prony_method: sorting abs of eigen-values done for imaginary parts'
 !write(6,*)
 !
 !!!!!!!!!!!! check the accurancy of Takagi factorization for imaginary
 !Hankel matirx !!!!!!!!!!!!
 write(6,*) 'prony_method: checking accuracy of Takagi factorization'
 allocate(cmattmp1(dimension_Hankel,dimension_Hankel),&
          cmattmp2(dimension_Hankel,dimension_Hankel), STAT=istat)
 allocate(cmattmp3(dimension_Hankel,dimension_Hankel),&
          cmattmp4(dimension_Hankel,dimension_Hankel), STAT=istat)
 cmattmp1 = czero
 cmattmp2 = czero
 cmattmp3 = czero
 cmattmp4 = czero
 error = 0.0d0
 !
 ! cmattmp1 = diagional matirx (hankel_imag_eigval)
 do ni = 1, dimension_Hankel
   cmattmp1(ni,ni) = zhankel_imag_eigval(ni)
   do nj = 1, dimension_Hankel
     cmattmp4(ni,nj) = zhankel_imag_eigvec(nj,ni)
   end do
 end do
 
 ! zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 ! C = alpha * A * B + beta * C
 ! cmattmp2 =  zhankel_imag_eigvec * hankel_imag_eigval 
 call zgemm('N', 'T', dimension_Hankel, dimension_Hankel, &
            dimension_Hankel, cunity, zhankel_imag_eigvec, &
            dimension_Hankel, cmattmp1, dimension_Hankel, &
            czero, cmattmp2, dimension_Hankel)
 ! 
 ! cmattmp3 =  zhankel_imag_eigvec * hankel_imag_eigval *
 ! zhankel_imag_eigvec.T
 call zgemm('N', 'T', dimension_Hankel, dimension_Hankel,&
            dimension_Hankel, cunity, cmattmp2, &
            dimension_Hankel, zhankel_imag_eigvec, &
            dimension_Hankel, czero, &
            cmattmp3, dimension_Hankel)
 ! 
 dtmp1 = 0.d0
 do ni = 1, dimension_Hankel
   do nj = 1, dimension_Hankel
     dtmp1 = max(dtmp1, cdabs( dcmplx(hankel_imag(ni,nj), 0.0d0)&
           - cmattmp3(ni,nj) ))
     error = error + cdabs( dcmplx(hankel_imag(ni,nj), 0.0d0)&
           - cmattmp3(ni,nj) )
   end do
 end do
 
 write(6,*)'prony_method: ', &
           'error of Takagi facotization for imagnary part',&
            error 
 write(6,*)'prony_method: the maxium error element is', dtmp1
 !write(6,*)
 
 deallocate(cmattmp1,cmattmp2,cmattmp3,cmattmp4, STAT=istat)
 
 !!!!!!!!!!!! calculate the roots of polynomial !!!!!!!!!!!!
 write(6,*) 'prony_method: calculate the root of polynomial' 
 allocate(cpara(dimension_Hankel - 1),&
          cpara_mat(dimension_Hankel - 1,dimension_Hankel - 1),&
          STAT=istat)
 allocate(prony_roots_imag(dimension_Hankel - 1), STAT=istat)
 cpara = czero
 cpara_mat = czero
 prony_roots_imag = czero
 
 do ni = 1, dimension_Hankel - 1
   cpara(ni) = zhankel_imag_eigvec(ni, prony_imag + 1)/&
   (zhankel_imag_eigvec(dimension_Hankel, prony_imag + 1))
   !if (ni <= 10 .or. (dimension_Hankel - 1 - ni) <= 10) then
   !  write(6,*)'numerator', &
   !  zhankel_imag_eigvec(ni, prony_imag + 1)
   !  write(6,*)'denominator',&
   !  (zhankel_imag_eigvec(dimension_Hankel,prony_imag + 1))
   !end if
 end do
 !
 do ni = 1, dimension_Hankel - 1
   do nj = 1, dimension_Hankel - 1
     if (nj == (ni + 1)) then
       cpara_mat(ni,nj) = cunity
     end if
   end do
   if (ni == dimension_Hankel - 1) then
     do nj = 1, dimension_Hankel - 1
       cpara_mat(ni,nj) = -cunity * cpara(nj)
     end do
   end if
 end do
 !
 lwork = 5*(dimension_Hankel-1)
 allocate(cmattmp1(dimension_Hankel - 1, dimension_Hankel - 1),&
          cmattmp2(dimension_Hankel - 1, dimension_Hankel - 1),&
          STAT=istat)
 allocate(cdwork(lwork), dwork(2*(dimension_Hankel-1)), STAT=istat)
 cdwork = czero
 cmattmp1 = czero
 cmattmp2 = czero
 dwork = 0.0d0
 dtmp1 = 0.0d0
 
 !call zgeev(JOBVL, JOBVR, N, A( lda, N ), LDA, W( N ), VL( ldvl, N ),
 !LDVL, VR( ldvr, N ), LDVR, WORK( LWORK ), LWORK, RWORK( 2N ), INFO)
 call zgeev('N','N', dimension_Hankel - 1, cpara_mat,&
            dimension_Hankel - 1, prony_roots_imag, &
      cmattmp1, dimension_Hankel - 1, cmattmp2, dimension_Hankel - 1,&
      cdwork, 5*(dimension_Hankel-1), dwork, info)
 
 deallocate(cdwork, dwork, STAT=istat)
 deallocate(cmattmp1,cmattmp2, STAT=istat)
 deallocate(zhankel_imag_eigvec, zhankel_imag_eigval, STAT=istat)
 write(6,*) 'prony_method: calculate the root of polynomial done for imaginary part' 
 !write(6,*)
 
 !!!!!!!!!!!! get the number of nonzero roots for imagniary part
 !!!!!!!!!!!!!
 write(6,*) 'prony_method: picking out nonzero roots'
 allocate(cvectmp1(dimension_Hankel), STAT=istat)
 cvectmp1 = czero
 ncount = 0
 do ni = 1, dimension_Hankel - 1
   if (cdabs(prony_roots_imag(ni)) >= 10E-16 ) then
     ncount = ncount + 1
     cvectmp1(ncount) = prony_roots_imag(ni)
   end if
 end do
 
 if (ncount == (dimension_Hankel - 1)) then
   write(6,*) 'prony_method: there is no zero element'
   write(6,*) 'prony_method: number of nonzero element', ncount
   !write(6,*)
 else if (ncount .ne. (dimension_Hankel - 1)) then
   deallocate(prony_roots_imag)
   allocate(prony_roots_imag(ncount), STAT=istat)
   write(6,*) 'prony_method: number of nonzero element', ncount
   !write(6,*)
   do ni = 1, ncount
     prony_roots_imag(ni) = cvectmp1(ni)
   end do
 end if
 
 deallocate(cvectmp1, STAT=istat)
 write(6,*) 'prony_method: picking out nonzero roots done for imagniary part'
 !write(6,*)
     
 !!!!!!!!!!!! get the order of dabs(prony_roots_imag) in
 !dabs(prony_roots_imag) sorted in the ascending order !!!!!!!!!!!! 
    
 allocate(nindextmp(ncount), STAT=istat)
 allocate(dvectmp1(ncount), dvectmp2(ncount), STAT=istat)
 dvectmp1 = 0.0d0
 dvectmp2 = 0.0d0
 nindextmp = 0
 do ni = 1, ncount
   dvectmp1(ni) = cdabs(prony_roots_imag(ni))
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
 
 allocate(cvectmp1(prony_imag), STAT=istat)
 cvectmp1 = czero
 
 do ni = 1, prony_imag
   cvectmp1(ni) = prony_roots_imag(nindextmp(ni))
 end do
 
 deallocate(prony_roots_imag, STAT=istat)
 
 !!!!!!!!!!!! calculate the dissipation rates of imagniary part
 !!!!!!!!!!!!!
 write(6,*) 'prony_method: calculate the dissipation rates of imagniary part' 
 allocate(prony_roots_imag(prony_imag),&
          diss_rate_prony_imag(prony_imag),&
          STAT=istat)
 prony_roots_imag = cvectmp1
 diss_rate_prony_imag = czero
 
 cvectmp1 = czero
 do ni = 1, prony_imag
   ctmp1 = dcmplx(atan( dimag(prony_roots_imag(ni)),&
           dble(prony_roots_imag(ni))), 0.0d0)
   cvectmp1(ni) = -2.d0 * (dimension_Hankel - 1) /&
                  tt_for_prony(tt_sample_for_prony) &
             * ( dlog( cdabs(prony_roots_imag(ni)) ) + eye * ctmp1)
 end do
 !
 diss_rate_prony_imag = cvectmp1
 deallocate(cvectmp1, STAT=istat)
 write(6,*) 'prony_method: calculate dissipation rates done for imagniary part' 
 !write(6,*)
 
 !!!!!!!!!!!! calculate the coefficients of imagniary part !!!!!!!!!!!!
 ! calculate the coefficients by the least squares method
 write(6,*) 'prony_method: calculate the coefficients of imagniary part' 
 allocate(cmattmp1(tt_sample_for_prony, prony_imag), STAT=istat)
 allocate(cmattmp2(prony_imag, tt_sample_for_prony), STAT=istat)
 allocate(cmattmp3(prony_imag, prony_imag), STAT=istat)
 allocate(cmattmp4(prony_imag, prony_imag), STAT=istat)
 allocate(cvectmp1(tt_sample_for_prony), STAT=istat)
 allocate(nvectmp1(prony_imag), STAT=istat)
 cmattmp1 = czero
 cmattmp2 = czero
 cmattmp3 = czero
 cmattmp4 = czero
 cvectmp1 = czero
 nvectmp1 = 0
 
 do ni = 1, tt_sample_for_prony
   do nj = 1, prony_imag
       cmattmp1(ni,nj) = prony_roots_imag(nj)**(ni - 1)
       cmattmp2(nj,ni) = prony_roots_imag(nj)**(ni - 1)
   end do
 end do
 
 ! cmattmp3 = cmattmp2 * cmattmp1 = omega_imag.T * omega_imag
 call zgemm('N', 'N', prony_imag, prony_imag,&
            tt_sample_for_prony, cunity, cmattmp2, &
            prony_imag, cmattmp1, tt_sample_for_prony,&
            czero, cmattmp3, prony_imag)
 
 ! cvectmp1 = cmattmp2 * dimag(res_prony) = omega_imag.T *
 ! dimag(res_prony)
 do ni = 1, prony_imag
   do nj = 1, tt_sample_for_prony
       cvectmp1(ni) = cvectmp1(ni) +& 
       cmattmp2(ni,nj) * dimag(res_prony(nj))
   end do
 end do
 
 ! inverse cmattmp3
 lwork = 5 * prony_imag
 allocate(cvectmp2(lwork), STAT=istat)
 cvectmp2 = czero
 cmattmp4 = cmattmp3
 
 ! calculate the coefficients by the least squares method
 write(6,*) 'prony_method: calculate coeffs by least squares'
 call zgetrf(prony_imag, prony_imag, &
             cmattmp3, prony_imag, nvectmp1,info)
 call zgetri(prony_imag, cmattmp3, &
             prony_imag, nvectmp1, cvectmp2, lwork, info) 
 ! 
 deallocate(nvectmp1, cvectmp2, STAT=istat)
 !
 allocate(cmattmp5(prony_imag, prony_imag), STAT=istat)
 !
 cmattmp5 = czero
 call zgemm('N', 'N', prony_imag, prony_imag,&
            prony_imag, cunity, cmattmp4, &
            prony_imag, cmattmp3, prony_imag,&
            czero, cmattmp5, prony_imag)
 allocate(eta_prony_imag(prony_imag), STAT=istat)
 eta_prony_imag = czero
 
 ! cmattmp3 * cvectmp1
 do ni = 1, prony_imag
   do nj = 1, prony_imag
       eta_prony_imag(ni) = eta_prony_imag(ni)&
          + cmattmp3(ni,nj) * cvectmp1(nj)
   end do
 end do
 
 deallocate(cmattmp1, cmattmp2, cmattmp3,&
            cmattmp4, cmattmp5, cvectmp1, STAT=istat)
 write(6,*) 'prony_method: calculate the coefficients done for imagniary part' 
 write(6,*)
 
 do ni = 1, prony_imag
     dtmp1 = dble(eta_prony_imag(ni))
     dtmp2 = dimag(eta_prony_imag(ni))
     if (dabs(dtmp1) <= 1.d-6) then
         dtmp1 = 0.d0
     end if
     if (dabs(dtmp2) <= 1.d-6) then
         dtmp2 = 0.d0
     end if
     eta_prony_imag(ni) = dcmplx(dtmp1,dtmp2)
 
    dtmp1 = dble(diss_rate_prony_imag(ni))
    dtmp2 = dimag(diss_rate_prony_imag(ni))
    if (dabs(dtmp1) <= 1.d-6) then
        dtmp1 = 0.d0
     end if
     if (dabs(dtmp2) <= 1.d-6) then
         dtmp2 = 0.d0
     end if
     diss_rate_prony_imag(ni) = dcmplx(dtmp1,dtmp2)
 end do

 write(6,*)'prony: output cb of imaginiary part'
 do nk = 1, (prony_imag)
   write(6,*)nk,dble(eta_prony_imag(nk)),&
   dimag(eta_prony_imag(nk))
 end do
 write(6,*)
 !
 write(6,*)'prony: output gamma of imagninary part'
 do nk = 1, (prony_imag)
   write(6,*)nk, -dble(diss_rate_prony_imag(nk)),&
   -dimag(diss_rate_prony_imag(nk))
 end do
 call flush(6)

 ntmp1 = 0
 do ni = 1, prony_imag
   if (dreal(eta_prony_imag(ni)) .eq. 0.d0 .and.&
       dimag(eta_prony_imag(ni)) .eq. 0.d0) then
     write(6,'(" prony: find zero cb and cd, No. ",I3)') ni
     cycle
   end if
   ntmp1 = ntmp1 + 1
 end do
 !
 if (ntmp1 .eq. prony_imag) then
    write(6,*) 'prony: all cb and cd are non-zero'
    goto 114
 end if
 !
 ntmp2 = 0
 allocate(eta_prony_tmp(ntmp1), &
          diss_rate_prony_tmp(ntmp1), STAT=istat)
 eta_prony_tmp = 0.d0
 diss_rate_prony_tmp = 0.d0
 do ni = 1, prony_imag
   if (dreal(eta_prony_imag(ni)) .eq. 0.d0 .and.&
       dimag(eta_prony_imag(ni)) .eq. 0.d0) then
     write(6,'(" prony: find zero cb and cd, No. ",I3)') ni
     cycle
   end if
   !
   ntmp2 = ntmp2 + 1
   eta_prony_tmp(ntmp2) = eta_prony_imag(ni)
   diss_rate_prony_tmp(ntmp2) = diss_rate_prony_imag(ni) 
 end do 
 prony_imag = ntmp1
 prony_imag_tmp = prony_imag
 write(6,*) 'number of rational cb and cd', prony_imag
 !
 deallocate(eta_prony_imag,diss_rate_prony_imag,STAT=istat)
 allocate(eta_prony_imag(prony_imag), &
          diss_rate_prony_imag(prony_imag),STAT=istat)
 eta_prony_imag = eta_prony_tmp
 diss_rate_prony_imag = diss_rate_prony_tmp
 deallocate(eta_prony_tmp,diss_rate_prony_tmp,STAT=istat)

 write(6,*)'prony: output rational cb of imaginiary part'
 do nk = 1, (prony_imag)
   write(6,*)nk,dble(eta_prony_imag(nk)),&
   dimag(eta_prony_imag(nk))
 end do
 write(6,*)
 !
 write(6,*)'prony: output rational gamma of imagninary part'
 do nk = 1, (prony_imag)
   write(6,*)nk, -dble(diss_rate_prony_imag(nk)),&
   -dimag(diss_rate_prony_imag(nk))
 end do
 call flush(6)
 114 continue
 
 if (lonly_imag .and. .not. lonly_real) then
   goto 1750
 end if

 1500 continue 
 !!!!!!!!!!!! Takagi factorization for real Hankel matirx !!!!!!!!!!!!
 write(6,*) 'prony_method: performing Takagi factorization for real part' 
 dtmp1 = 0.d0
 dtmp2 = 0.d0
 dtmp3 = 0.d0
 do ni = 1, dimension_Hankel
   do nj = 1, dimension_Hankel
       dtmp1 = dabs(hankel_real(ni,nj) - hankel_real(nj,ni))
       dtmp2 = dtmp2 + dtmp1
       dtmp3 = max(dtmp1, dtmp3)
   end do
 end do
 
 write(6,*)'prony_method: sum over Hankel_real - Hankel_real.T', dtmp2
 write(6,*)'prony_method: max error of Hankel_real - Hankel_real.T', dtmp3
 !write(6,*)
 
 allocate(hankel_tmp_eigval(dimension_Hankel),&
          hankel_tmp_eigvec(dimension_Hankel,dimension_Hankel),&
          STAT=istat)
 hankel_tmp_eigvec = hankel_real
 hankel_tmp_eigval = 0.0d0
 ! 
 lwork = (1 + 6 * dimension_Hankel + 2 * dimension_Hankel**2) * 2
 liwork = (3 + 5 * dimension_Hankel) * 2
!!! DEBUG
 allocate(dwork(lwork), iwork(liwork), STAT=istat)    
!!! DEBUG
 dwork = 0.0d0  
 iwork = 0
 
 !call dsyevd(JOBZ, UPLO, N, A( lda, N ), LDA, W( N ), WORK( LWORK ),
 !LWORK, IWORK, LIWORK, INFO)
 write(6,*) 'prony_method: solving eigen-value of Hankel matrix'
 call dsyevd('V','L', dimension_Hankel, hankel_tmp_eigvec,&
             dimension_Hankel, hankel_tmp_eigval,&
             dwork, lwork, iwork, liwork, info)
 write(6,*) 'prony_method: eigen-value calculation done'
 if (info .ne. 0) then
  write(6,*)'prony_method:',& 
           ' error! eigenvalue calculation', &
           ' failed for hankel_real', info
  stop
 end if
 
 deallocate(dwork, iwork, STAT=istat)
 
 allocate(nindextmp(dimension_Hankel), STAT=istat)
 allocate(dvectmp1(dimension_Hankel), STAT=istat)
 allocate(cmattmp1(dimension_Hankel,dimension_Hankel),&
          cmattmp2(dimension_Hankel,dimension_Hankel), STAT=istat)
 allocate(zhankel_real_eigval(dimension_Hankel),&
          zhankel_real_eigvec(dimension_Hankel,dimension_Hankel),&
          STAT=istat)
!
 nindextmp = 0
 dvectmp1 = 0.0d0
 cmattmp1 = czero
 cmattmp2 = czero
 
 do nj = 1, dimension_Hankel
   if (hankel_tmp_eigval(nj) < 0.0d0) then
     do ni = 1, dimension_Hankel
       zhankel_real_eigvec(ni,nj) = -eye *&
        hankel_tmp_eigvec(ni,nj)
     end do
   else if (hankel_tmp_eigval(nj) > 0.0d0) then
     do ni = 1, dimension_Hankel
       zhankel_real_eigvec(ni,nj) = cunity *&
        hankel_tmp_eigvec(ni,nj)
     end do
   end if
 end do
 
 do ni = 1, dimension_Hankel
   dvectmp1(ni) = dabs(hankel_tmp_eigval(ni))
   nindextmp(ni) = ni
 end do
 
 deallocate(hankel_tmp_eigval, hankel_tmp_eigvec, STAT=istat)
 write(6,*) 'prony_method: Takagi factorization done for real part' 
 !write(6,*)
 
 cmattmp1 = zhankel_real_eigvec
 
 ! sort dabs(hankel_real_eigval) in the descending order
 write(6,*) 'prony_method: sorting abs of eigen-values'
 do ni = 1, dimension_Hankel -1
   do nj = ni, 1, -1
     if (dvectmp1(nj) < dvectmp1(nj + 1)) then
       dtmp1 = dvectmp1(nj + 1)
       dvectmp1(nj + 1) = dvectmp1(nj)
       dvectmp1(nj) = dtmp1
       
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
 
 zhankel_real_eigval = dcmplx(dvectmp1, 0.0d0)
 zhankel_real_eigvec = cmattmp2
 
 deallocate(nindextmp, STAT=istat)
 deallocate(dvectmp1, STAT=istat)
 deallocate(cmattmp1, cmattmp2, STAT=istat)
 write(6,*) 'prony_method: sorting abs of eigen-values done for imaginary parts'
 !write(6,*)
 
 !!!!!!!!!!!! check the accurancy of Takagi factorization for real
 !Hankel matirx !!!!!!!!!!!!
 write(6,*) 'prony_method: checking accuracy of Takagi factorization'
 allocate(cmattmp1(dimension_Hankel,dimension_Hankel),&
          cmattmp2(dimension_Hankel,dimension_Hankel), STAT=istat)
 allocate(cmattmp3(dimension_Hankel,dimension_Hankel),&
          cmattmp4(dimension_Hankel,dimension_Hankel), STAT=istat)
 cmattmp1 = czero
 cmattmp2 = czero
 cmattmp3 = czero
 cmattmp4 = czero
 error = 0.0d0
 
 ! cmattmp1 = diagional matirx (hankel_real_eigval)
 do ni = 1, dimension_Hankel
   cmattmp1(ni,ni) = zhankel_real_eigval(ni)
   do nj = 1, dimension_Hankel
     cmattmp4(ni,nj) = zhankel_real_eigvec(nj,ni)
   end do
 end do
 
 ! zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 ! C = alpha * A * B + beta * C
 ! cmattmp2 =  zhankel_real_eigvec * hankel_real_eigval 
 call zgemm('N', 'T', dimension_Hankel, dimension_Hankel,&
            dimension_Hankel, cunity, zhankel_real_eigvec,&
            dimension_Hankel, cmattmp1, dimension_Hankel,&
            czero, cmattmp2, dimension_Hankel)
  
 ! cmattmp3 =  zhankel_real_eigvec * hankel_real_eigval *
 ! zhankel_real_eigvec.T
 call zgemm('N', 'T', dimension_Hankel, dimension_Hankel,&
            dimension_Hankel, cunity, cmattmp2, &
            dimension_Hankel, zhankel_real_eigvec,& 
            dimension_Hankel, czero, cmattmp3, dimension_Hankel)
 dtmp1 = 0.0d0
 do ni = 1, dimension_Hankel
   do nj = 1, dimension_Hankel
     dtmp1 = max(dtmp1, cdabs( dcmplx(hankel_real(ni,nj),&
     0.0d0) - cmattmp3(ni,nj) ))
     error = error + cdabs( dcmplx(hankel_real(ni,nj),&
     0.0d0) - cmattmp3(ni,nj) )
   end do
 end do
 
 write(6,*)'prony_method:',& 
          ' error of Takagi facotization for real part',&
            error 
 write(6,*)'prony_method: the maxium error element is', dtmp1
 !write(6,*)
 
 deallocate(cmattmp1,cmattmp2,cmattmp3,cmattmp4, STAT=istat)
 
 !!!!!!!!!!!! calculate the roots of polynomial !!!!!!!!!!!!
 write(6,*) 'prony_method: calculate the root of polynomial' 
 allocate(cpara(dimension_Hankel - 1),&
          cpara_mat(dimension_Hankel - 1,&
          dimension_Hankel - 1), STAT=istat)
 allocate(prony_roots_real(dimension_Hankel - 1), STAT=istat)
 cpara = czero
 cpara_mat = czero
 prony_roots_real = czero
 
 do ni = 1, dimension_Hankel - 1
   cpara(ni) = zhankel_real_eigvec(ni, prony_real&
    + 1)/(zhankel_real_eigvec(dimension_Hankel, prony_real + 1))
   !if (ni <= 10) then
   !  write(6,*)'numerator', zhankel_real_eigvec(ni, prony_real + 1)
   !  write(6,*)'denominator',&
   !         (zhankel_real_eigvec(dimension_Hankel, prony_real + 1))
   !end if
 end do
 
 ntmp1 = 0
 do ni = 1, dimension_Hankel - 1
   if (isnan(cdabs(cpara(ni)))) then
    ntmp1 = ntmp1 + 1 
   end if
 end do
 write(6,*)'prony_method: number of nan in cpara', ntmp1
 
 do ni = 1, dimension_Hankel - 1
   do nj = 1, dimension_Hankel - 1
     if (nj == (ni + 1)) then
       cpara_mat(ni,nj) = cunity
     end if
   end do
   if (ni == dimension_Hankel - 1) then
     do nj = 1, dimension_Hankel - 1
       cpara_mat(ni,nj) = -cunity * cpara(nj)
     end do
   end if
 end do
 
 lwork = 5*(dimension_Hankel-1)
 allocate(cmattmp1(dimension_Hankel - 1, & 
                   dimension_Hankel - 1),&
          cmattmp2(dimension_Hankel - 1, & 
                   dimension_Hankel - 1), STAT=istat)
 allocate(cdwork(lwork), dwork(2*(dimension_Hankel-1)), STAT=istat)
 cdwork = czero
 cmattmp1 = czero
 cmattmp2 = czero
 dwork = 0.0d0
 dtmp1 = 0.0d0
 
 
 !call zgeev(JOBVL, JOBVR, N, A( lda, N ), LDA, W( N ), VL( ldvl, N ),
 !LDVL, VR( ldvr, N ), LDVR, WORK( LWORK ), LWORK, RWORK( 2N ), INFO)
 call zgeev('N','N', dimension_Hankel - 1, cpara_mat,&
            dimension_Hankel - 1, prony_roots_real, &
            cmattmp1, dimension_Hankel - 1, cmattmp2,&
            dimension_Hankel - 1, cdwork,&
            5*(dimension_Hankel-1), dwork, info)
 
 deallocate(cdwork, dwork, STAT=istat)
 deallocate(cmattmp1,cmattmp2, STAT=istat)
 deallocate(zhankel_real_eigvec, zhankel_real_eigval, STAT=istat)
 write(6,*) 'prony_method: calculate the root of polynomial done for real part' 
 !write(6,*)
 
 !!!!!!!!!!!! get the number of nonzero roots for real part !!!!!!!!!!!!
 
 allocate(cvectmp1(dimension_Hankel), STAT=istat)
 cvectmp1 = czero
 ncount = 0
 do ni = 1, dimension_Hankel - 1
   if (cdabs(prony_roots_real(ni)) >= 10E-16 ) then
     ncount = ncount + 1
     cvectmp1(ncount) = prony_roots_real(ni)
   end if
 end do
 
 if (ncount == (dimension_Hankel - 1)) then
   write(6,*) 'prony_method: there is no zero element for real part'
   write(6,*) 'prony_method: number of nonzero element', ncount
   write(6,*)
 else if (ncount .ne. (dimension_Hankel - 1)) then
   deallocate(prony_roots_real)
   allocate(prony_roots_real(ncount), STAT=istat)
   write(6,*) 'prony_method:',&
         ' number of nonzero element for real part', ncount
   write(6,*)
   do ni = 1, ncount
     prony_roots_real(ni) = cvectmp1(ni)
   end do
 end if
 
 deallocate(cvectmp1, STAT=istat)
     
 !!!!!!!!!!!! get the order of dabs(prony_roots_real) in
 !dabs(prony_roots_real) sorted in the ascending order !!!!!!!!!!!! 
    
 allocate(nindextmp(ncount), STAT=istat)
 allocate(dvectmp1(ncount), dvectmp2(ncount), STAT=istat)
 dvectmp1 = 0.0d0
 dvectmp2 = 0.0d0
 nindextmp = 0
 do ni = 1, ncount
   dvectmp1(ni) = cdabs(prony_roots_real(ni))
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
       
       ntmp1 = nindextmp(nj)
       nindextmp(nj) = nindextmp(nj + 1)
       nindextmp(nj + 1) = ntmp1
     end if
   end do
 end do

 
 allocate(cvectmp1(prony_real), STAT=istat)
 cvectmp1 = czero
 
 do ni = 1, prony_real
   cvectmp1(ni) = prony_roots_real(nindextmp(ni))
 end do
 
 deallocate(prony_roots_real, STAT=istat)
 
 !!!!!!!!!!!! calculate the dissipation rates of real part !!!!!!!!!!!!
 write(6,*) 'prony_method: calculate the dissipation rates of real part' 
 allocate(prony_roots_real(prony_real),&
          diss_rate_prony_real(prony_real), STAT=istat)
 prony_roots_real = cvectmp1
 diss_rate_prony_real = czero
 
 cvectmp1 = czero
 do ni = 1, prony_real
   ctmp1 = dcmplx(atan( dimag(prony_roots_real(ni)),&
           dble(prony_roots_real(ni))), 0.0d0)
   cvectmp1(ni) = -2.d0 * (dimension_Hankel - 1) /&
               tt_for_prony(tt_sample_for_prony) &
             * ( dlog( cdabs(prony_roots_real(ni)) ) + eye * ctmp1)
 end do
 
 diss_rate_prony_real = cvectmp1
 deallocate(cvectmp1, STAT=istat)
 write(6,*) 'prony_method: calculate dissipation rates done for real part' 
 !write(6,*)
 
 !!!!!!!!!!!! calculate the coefficients of real part !!!!!!!!!!!!
 ! calculate the coefficients by the least squares method
 write(6,*) 'prony_method: calculate the coefficients of real part' 
 allocate(cmattmp1(tt_sample_for_prony, prony_real), STAT=istat)
 allocate(cmattmp2(prony_real, tt_sample_for_prony), STAT=istat)
 allocate(cmattmp3(prony_real, prony_real), STAT=istat)
 allocate(cmattmp4(prony_real, prony_real), STAT=istat)
 allocate(cvectmp1(tt_sample_for_prony), STAT=istat)
 allocate(nvectmp1(prony_real), STAT=istat)
 cmattmp1 = czero
 cmattmp2 = czero
 cmattmp3 = czero
 cmattmp4 = czero
 cvectmp1 = czero
 nvectmp1 = 0
 
 do ni = 1, tt_sample_for_prony
   do nj = 1, prony_real
     cmattmp1(ni,nj) = prony_roots_real(nj)**(ni - 1)
     cmattmp2(nj,ni) = prony_roots_real(nj)**(ni - 1)
   end do
 end do
 
 ! cmattmp3 = cmattmp2 * cmattmp1 = omega_real.T * omega_real
 call zgemm('N', 'N', prony_real, prony_real,&
            tt_sample_for_prony, cunity, cmattmp2, &
            prony_real, cmattmp1, tt_sample_for_prony,&
            czero, cmattmp3, prony_real)
 
 ! cvectmp1 = cmattmp2 * dble(res_prony) = omega_real.T *
 ! dble(res_prony)
 do ni = 1, prony_real
   do nj = 1, tt_sample_for_prony
       cvectmp1(ni) = cvectmp1(ni) &
     + cmattmp2(ni,nj) * dble(res_prony(nj))
   end do
 end do
 
 ! inverse cmattmp3
 lwork = 5 * prony_real
 allocate(cvectmp2(lwork), STAT=istat)
 cvectmp2 = czero
 cmattmp4 = cmattmp3
 
 ! calculate the coefficients by the least squares method
 call zgetrf(prony_real, prony_real, cmattmp3,&
             prony_real, nvectmp1, info)
 call zgetri(prony_real, cmattmp3, prony_real,&
             nvectmp1, cvectmp2, lwork, info) 
 
 deallocate(nvectmp1, cvectmp2, STAT=istat)
 
 allocate(cmattmp5(prony_real, prony_real), STAT=istat)
 
 cmattmp5 = czero
 call zgemm('N', 'N', prony_real, prony_real, prony_real,&
            cunity, cmattmp4, prony_real, cmattmp3,&
            prony_real, czero, cmattmp5, prony_real)
 
 allocate(eta_prony_real(prony_real), STAT=istat)
 eta_prony_real = czero
 
 ! cmattmp3 * cvectmp1
 do ni = 1, prony_real
   do nj = 1, prony_real
     eta_prony_real(ni) = eta_prony_real(ni)&
                   + cmattmp3(ni,nj) * cvectmp1(nj)
   end do
 end do
 
 deallocate(cmattmp1, cmattmp2, cmattmp3,&
            cmattmp4, cmattmp5, cvectmp1, STAT=istat)
 write(6,*) 'prony_method: calculate the coefficients done for real part' 
 write(6,*)
 
 do ni = 1, prony_real
     dtmp1 = dble(eta_prony_real(ni))
     dtmp2 = dimag(eta_prony_real(ni))
     if (dabs(dtmp1) <= 1.d-6) then
         dtmp1 = 0.d0
     end if
     if (dabs(dtmp2) <= 1.d-6) then
         dtmp2 = 0.d0
     end if
     eta_prony_real(ni) = dcmplx(dtmp1,dtmp2)
 
     dtmp1 = dble(diss_rate_prony_real(ni))
     dtmp2 = dimag(diss_rate_prony_real(ni))
     if (dabs(dtmp1) <= 1.d-6) then
         dtmp1 = 0.d0
     end if
     if (dabs(dtmp2) <= 1.d-6) then
         dtmp2 = 0.d0
     end if
     diss_rate_prony_real(ni) = dcmplx(dtmp1,dtmp2)
 end do

 write(6,*)'prony: output cb of real part'
 do nk = 1, (prony_real)
   write(6,*)nk,dble(eta_prony_real(nk)),&
   dimag(eta_prony_real(nk))
 end do
 write(6,*)
 !
 write(6,*)'prony: output gamma of real part'
 do nk = 1, (prony_real)
   write(6,*)nk, -dble(diss_rate_prony_real(nk)),&
   -dimag(diss_rate_prony_real(nk))
 end do
 call flush(6)

 ntmp1 = 0
 do ni = 1, prony_real
   if (dreal(eta_prony_real(ni)) .eq. 0.d0 .and.&
       dimag(eta_prony_real(ni)) .eq. 0.d0) then
     write(6,'(" prony: find zero cb and cd, No. ",I3)') ni
     cycle
   end if
   ntmp1 = ntmp1 + 1
 end do
 !
 if (ntmp1 .eq. prony_real) then
    write(6,*) 'prony: all cb and cd are non-zero'
    goto 514
 end if
 !
 ntmp2 = 0
 allocate(eta_prony_tmp(ntmp1), &
          diss_rate_prony_tmp(ntmp1), STAT=istat)
 eta_prony_tmp = 0.d0
 diss_rate_prony_tmp = 0.d0
 do ni = 1, prony_real
   if (dreal(eta_prony_real(ni)) .eq. 0.d0 .and.&
       dimag(eta_prony_real(ni)) .eq. 0.d0) then
     write(6,'(" prony: find zero cb and cd, No. ",I3)') ni
     cycle
   end if
   !
   ntmp2 = ntmp2 + 1
   eta_prony_tmp(ntmp2) = eta_prony_real(ni)
   diss_rate_prony_tmp(ntmp2) = diss_rate_prony_real(ni) 
 end do 
 prony_real = ntmp1
 prony_real_tmp = prony_real
 write(6,*) 'number of rational cb and cd', prony_real
 !
 deallocate(eta_prony_real,diss_rate_prony_real,STAT=istat)
 allocate(eta_prony_real(prony_real), &
          diss_rate_prony_real(prony_real),STAT=istat)
 eta_prony_real = eta_prony_tmp
 diss_rate_prony_real = diss_rate_prony_tmp
 deallocate(eta_prony_tmp,diss_rate_prony_tmp,STAT=istat)

 write(6,*)'prony: output rational cb of real part'
 do nk = 1, (prony_real)
   write(6,*)nk,dble(eta_prony_real(nk)),&
   dimag(eta_prony_real(nk))
 end do
 write(6,*)
 !
 write(6,*)'prony: output rational gamma of real part'
 do nk = 1, (prony_real)
   write(6,*)nk, -dble(diss_rate_prony_real(nk)),&
   -dimag(diss_rate_prony_real(nk))
 end do
 call flush(6)
 514 continue
 
 !!!!!!!!!!!! output dissipation rate and coefficients !!!!!!!!!!!!
 1750 continue
 write(6,*) 'prony_method: output dissipation rate and coefficients' 
 allocate(eta_prony(prony_real + prony_imag),&
          diss_rate_prony(prony_real + prony_imag), STAT=istat)
 do ni = 1, (prony_real + prony_imag)
   if (ni <= prony_real) then
     eta_prony(ni) = eta_prony_real(ni)
     diss_rate_prony(ni) = diss_rate_prony_real(ni)
   else if (ni > prony_real) then
     eta_prony(ni) = eta_prony_imag(ni - prony_real) * eye
     diss_rate_prony(ni) = diss_rate_prony_imag(ni - prony_real)
   end if
 end do
 
 !do ni = 1, (prony_real + prony_imag)
 !  dtmp1 = dble(eta_prony(ni))
 !  dtmp2 = dimag(eta_prony(ni))
 !  if (dabs(dtmp1) <= 1.d-8) then
 !    dtmp1 = 0.d0
 !  end if
 !  if (dabs(dtmp2) <= 1.d-8) then
 !    dtmp2 = 0.d0
 !  end if
 !  eta_prony(ni) = dcmplx(dtmp1,dtmp2)
 !  !
 !  dtmp1 = dble(diss_rate_prony(ni))
 !  dtmp2 = dimag(diss_rate_prony(ni))
 !  if (dabs(dtmp1) <= 1.d-8) then
 !    dtmp1 = 0.d0
 !  end if
 !  if (dabs(dtmp2) <= 1.d-8) then
 !    dtmp2 = 0.d0
 !  end if
 !  diss_rate_prony(ni) = dcmplx(dtmp1,dtmp2)
 !end do
 
 !open(unit = 90, file='eta_real.data')
 !rewind(90)
 !open(unit = 91, file='eta_imag.data')
 !rewind(91)
 !open(unit = 92, file='rate_real.data')
 !rewind(92)
 !open(unit = 93, file='rate_imag.data')
 !rewind(93)
 !
 !do ni =1, (prony_real + prony_imag)
 !    write(90,*) dble(eta_prony(ni))
 !    write(91,*) dimag(eta_prony(ni))
 !    write(92,*) dble(diss_rate_prony(ni))
 !    write(93,*) dimag(diss_rate_prony(ni))
 !end do
 !
 !close(90)
 !close(91)
 !close(92)
 !close(93)
 
 ! only valid for bcenter = 0
 open(unit = 94, file='PRONY_CB.data')
 rewind(94)
 write(94,*)'icor, Re[cb], Imag[cb]'
 open(unit = 95, file='PRONY_CD.data')
 rewind(95)
 write(95,*)'icor, Re[cd], Imag[cd]'
 open(unit = 96, file='PRONY_GAMA.data')
 rewind(96)
 write(96,*)'icor, Re[cgama], Imag[cgama]'
 do nk = 1, (prony_real + prony_imag)
   write(94,1001)nk,dble(eta_prony(nk))/hbar**2,&
    dimag(eta_prony(nk))/hbar**2
   write(95,1001)nk, dble(eta_prony(nk))/hbar**2,&
    -dimag(eta_prony(nk))/hbar**2
   write(96,1001)nk, -dble(diss_rate_prony(nk))/hbar,&
    -dimag(diss_rate_prony(nk))/hbar
 end do
 close(94)
 close(95)
 close(96)
  
 write(6,*)'prony_method: output cb'
 do nk = 1, (prony_real + prony_imag)
   write(6,1001)nk,dble(eta_prony(nk)),&
    dimag(eta_prony(nk))
 end do
 !write(6,*)'prony_method: out put cb'
 !write(6,*) (dreal(eta_prony(ni)),&
 !            dimag(eta_prony(ni)),&
 !            ni = 1,(prony_real + prony_imag))
 write(6,*)    
 write(6,*)'prony_method: output gamma'
 do nk = 1, (prony_real + prony_imag)
   write(6,1001)nk, -dble(diss_rate_prony(nk)),&
    -dimag(diss_rate_prony(nk))
 end do
 !write(6,*) (dreal(diss_rate_prony(ni)),&
 !            dimag(diss_rate_prony(ni)),&
 !            ni = 1,(prony_real + prony_imag))
 write(6,*)
 
 1001 format(1(I3,1x),2(1x, f37.32))
 1000 format(20(f17.13,2x))

 2000 continue
end subroutine prony_method
