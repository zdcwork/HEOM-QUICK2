!subroutine prony_job(sys_bath_coupling_tmp, band_width_tmp, band_center_tmp, T_tmp, sigma_tmp, npade_tmp, nprony_real_tmp, nprony_imag_tmp, &
!dimension_Hankel_tmp, tt_max_for_prony_tmp, output_eta, output_diss_rate)
program prony_job
! input  sys_bath_coupling, band_width, band_center, T and sigma
! output cb and diss_rate calculated by prony

    implicit none
    
    !real*8, intent (in)      :: sys_bath_coupling_tmp, band_width_tmp, band_center_tmp, T_tmp
    !integer, intent (in)     :: sigma_tmp
    !integer, intent (in)     :: npade_tmp, nprony_real, nprony_imag, dimension_Hankel_tmp, tt_max_for_prony_tmp
    !complex*16, intent (out)  :: output_eta(nprony_tmp), output_diss_rate(nprony_tmp)
    
    integer    :: ni, nj, nk, nl, nm, ncount, nsign
    integer    :: npade 
    integer    :: dimension_Hankel, tt_sample_for_prony
    integer    :: ntmp1, ntmp2, ntmp3
    integer, allocatable    :: nvectmp1(:)
    integer, allocatable    :: nindextmp(:)
    
    real*8     :: hbar, pi
    real*8     :: sys_bath_coupling, band_width, band_center, chemical_potential, T, sigma
    real*8     :: tt_min_for_prony, tt_step_for_prony
    real*8     :: tt_max_for_prony
    real*8     :: dtmp1, dtmp2, dtmp3, error
    real*8, allocatable    :: tt_for_prony(:)
    real*8, allocatable :: dvectmp1(:), dvectmp2(:)
    
    
    integer  :: nprony_real, nprony_imag, istat, info, lwork, liwork
    integer, allocatable  :: iwork(:)
    complex*16 :: cunit, eye, czero
    complex*16 :: ctmp1, ctmp2, ctmp3
    real*8, allocatable    :: hankel_real(:,:), hankel_imag(:,:)
    real*8, allocatable    :: hankel_tmp_eigval(:), hankel_tmp_eigvec(:,:)
    real*8, allocatable    :: dwork(:)
    complex*16, allocatable :: zpeigv(:), zpcoef(:), zppole_tmp(:,:), zppole(:), cpara(:), cpara_mat(:,:)
    complex*16, allocatable :: eta_psd(:), eta_psd_dag(:), diss_rate_psd(:)
    complex*16, allocatable :: eta_prony(:), diss_rate_prony(:), res_prony(:)
    complex*16, allocatable :: cvectmp1(:), cvectmp2(:), cdwork(:)
    complex*16, allocatable :: cmattmp1(:,:), cmattmp2(:,:), cmattmp3(:,:), cmattmp4(:,:), cmattmp5(:,:)
    complex*16, allocatable :: zhankel_real_eigval(:), zhankel_imag_eigval(:), zhankel_real_eigvec(:,:), zhankel_imag_eigvec(:,:)
    complex*16, allocatable :: prony_roots_real(:), prony_roots_imag(:)
    complex*16, allocatable :: diss_rate_prony_real(:), diss_rate_prony_imag(:)
    complex*16, allocatable :: eta_prony_real(:), eta_prony_imag(:)
    
    !real*8, external        :: datan
    
    !!!!!!!!!!!! necessary and important constant !!!!!!!!!!!!
    
    hbar = 0.658211928d0
    pi = 3.14159265358979323846d0
    cunit = (1.0d0, 0.0d0)
    eye = (0.0d0, 1.0d0)
    czero = (0.0d0, 0.0d0)
    
    nsign = 2
    tt_min_for_prony = 0.0d0
    
    !!!!!!!!!!!! input parameters !!!!!!!!!!!!
    
    npade = 650
    sigma = 1.d0
    
    nprony_real = 1
    nprony_imag = 11
    dimension_Hankel = 3250      !!! default value = 1500
    tt_max_for_prony = 184.0d0    !!! tt_max_for_prony should be sufficiently large to describe long time dynamics
    
    sys_bath_coupling = 1.0d0       !eV, is equal to linewidth(norbs,nalf) in HEOM-QUICK
    band_width = 5.0d0           !eV, is equal to gammap(nvar2,nspin,malf) or bandwidth(nalf) in HEOM-QUICK
    band_center = 0.0d0          !eV, is equal to dwp(iorbs2,ispin,ialf) or bcenter(nalf,nspin) in HEOM-QUICK
    T = 0.001d0               !eV, is equal to dinvbeta(nalf) in HEOM-QUICK
    !chemical_potential = 0.0d0
    
    !sigma = dble(sigma_tmp)
    !npade = int(npade_tmp)
    
    !nprony_real = nprony_real_tmp
    !nprony_imag = nprony_imag_tmp
    !dimension_Hankel = int(dimension_Hankel_tmp)      !!! default value = 3250
    !tt_max_for_prony = dble(tt_max_for_prony_tmp)     !!! default value = 180.0d0 Notice that tt_max_for_prony should not be too large
    
    !sys_bath_coupling = dble(sys_bath_coupling_tmp)   !eV, is equal to linewidth(norbs,nalf) in HEOM-QUICK
    !band_width = dble(band_width_tmp)             !eV, is equal to gammap(nvar2,nspin,malf) or bandwidth(nalf) in HEOM-QUICK
    !band_center = dble(band_center_tmp)            !eV, is equal to dwp(nvar2,nspin,nalf) or bcenter(nalf,nspin) in HEOM-QUICK
    !T = dble(T_tmp)                         !eV, is equal to dinvbeta(nalf) in HEOM-QUICK
    
    !!!!!!!!!!!! get pade poles and coefficients !!!!!!!!!!!!
    
    allocate(zpeigv(npade), zpcoef(npade), zppole(npade), zppole_tmp(npade,nsign), STAT=istat)
    allocate(eta_psd(npade + 1), eta_psd_dag(npade + 1), diss_rate_psd(npade + 1), STAT=istat)
    zpeigv = czero
    zpcoef = czero
    zppole = czero
    zppole_tmp = czero
    
    call zpout_psd(npade, nsign, zpeigv, zpcoef, zppole_tmp)
    
    if (sigma == 1.0) then
        do ni = 1, npade
            zppole(ni) = zppole_tmp(ni,1)
        end do
    else if (sigma == -1.0) then
        do ni = 1, npade
            zppole(ni) = zppole_tmp(ni,2)
        end do
    else
        write(6,*)'prony_job: wrong sigma!', sigma
        stop
    end if
    
    deallocate(zppole_tmp, STAT=istat)
    
    diss_rate_psd(1) = (band_width - eye * sigma * band_center) / hbar
    ctmp1 = eye * diss_rate_psd(1) / T * hbar
    ctmp2 = dcmplx(5.d-1, 0.d0)
    do ni=1, npade
        ctmp2 = ctmp2 + zpcoef(ni) / (ctmp1 + zppole(ni)) + &
                   zpcoef(ni) / (ctmp1 - zppole(ni))
    end do
    eta_psd(1) = 0.5d0 * sys_bath_coupling * band_width * ctmp2 / hbar**2

    eta_psd_dag(1) = dconjg(eta_psd(1))
    do ni = 2, npade + 1
         diss_rate_psd(ni) = ( -1.d0 * eye * sigma * T * zppole(ni-1) ) / hbar
         eta_psd(ni) = 2.d0 * eye * zpcoef(ni-1) * sys_bath_coupling * 0.5d0 * band_width**2 * T &
                       / ( (zppole(ni-1) * T - band_center)**2 + band_width**2 ) / hbar**2
         eta_psd_dag(ni) = dconjg(eta_psd(ni))
    end do
    
    write(6,*) 'diss_rate_psd for simga = 1'
    write(6,*) (-1.d0 * diss_rate_psd(ni) * hbar, ni = 1, 10)
    write(6,*) 'eta_psd for simga = 1'
    write(6,*) (eta_psd(ni), ni = 1, 10)
    write(6,*)

    deallocate(zpeigv, zpcoef, zppole, STAT=istat)
    
    
    !!!!!!!!!!!! build real and imaginary Hankel matirx for C(t) !!!!!!!!!!!!
    
    dimension_Hankel = dimension_Hankel + 1
    tt_sample_for_prony = 2 * dimension_Hankel - 1
    
    
    if (tt_max_for_prony < tt_min_for_prony) then
        dtmp1 = tt_max_for_prony
        tt_max_for_prony = tt_min_for_prony
        tt_min_for_prony = dtmp1
    else if (tt_max_for_prony == tt_min_for_prony) then
        write(6,*) 'wrong initial and final time', tt_min_for_prony, tt_max_for_prony
    end if
    
    tt_step_for_prony = (tt_max_for_prony - tt_min_for_prony) / (tt_sample_for_prony - 1)
    allocate(tt_for_prony(tt_sample_for_prony), res_prony(tt_sample_for_prony), STAT=istat)
    res_prony = czero
    
    do ni = 1, tt_sample_for_prony
        tt_for_prony(ni) = tt_min_for_prony + dble(ni - 1) * tt_step_for_prony
    end do
    
    write(6,*) 'prony_job: initial time =', tt_for_prony(1)
    write(6,*) 'prony_job: time step   =', tt_for_prony(2)
    write(6,*) 'prony_job: final time  =', tt_for_prony(tt_sample_for_prony)
    write(6,*)
    
    !if (tt_step_for_prony <= 0.06d0 .and. tt_step_for_prony >= 0.005d0) then
    !    write(6,*) 'prony_job: time step is ok'
    !    write(6,*)
    !else if (tt_step_for_prony <= 0.005d0 .and. tt_step_for_prony >= 0.0d0) then
    !    write(6,*) 'prony_job: time step is too small'
    !    write(6,*)
    !else if (tt_step_for_prony > 0.06d0) then
    !    write(6,*) 'prony_job: time step is too large'
    !    write(6,*)
    !else if (tt_step_for_prony <= 0.0d0 .or. tt_max_for_prony <= 0.0d0) then
    !    write(6,*) 'prony_job: wrong time information!', tt_step_for_prony, tt_max_for_prony
    !    stop
    !end if
    
    do ni = 1, tt_sample_for_prony
        do nj = 1, npade + 1
            res_prony(ni) = res_prony(ni) + eta_psd(nj) * cdexp(-1.d0 * diss_rate_psd(nj) * tt_for_prony(ni))
        end do
    end do
    
    allocate(hankel_real(dimension_Hankel, dimension_Hankel), hankel_imag(dimension_Hankel, dimension_Hankel), STAT=istat)
    hankel_real = 0.d0
    hankel_imag = 0.d0
    do ni = 1, dimension_Hankel
        do nj = 1, dimension_Hankel
            hankel_real(ni,nj) = dble(res_prony(ni + nj - 1))
            hankel_imag(ni,nj) = dimag(res_prony(ni + nj - 1))
        end do
    end do
     
    ! add the symmetry check at here!
    
    !!!!!!!!!!!! Takagi factorization for imaginary Hankel matirx !!!!!!!!!!!!
    
    allocate(hankel_tmp_eigval(dimension_Hankel), hankel_tmp_eigvec(dimension_Hankel,dimension_Hankel), STAT=istat)
    hankel_tmp_eigvec = hankel_imag
    hankel_tmp_eigval = 0.d0
    
    lwork = (1 + 6 * dimension_Hankel + 2 * dimension_Hankel**2) * 2
    liwork = (3 + 5 * dimension_Hankel) * 2
    allocate(dwork(lwork), iwork(liwork), STAT=istat)
    dwork = 0.0d0
    iwork = 0
    !call dsyevd(JOBZ, UPLO, N, A( lda, N ), LDA, W( N ), WORK( LWORK ), LWORK, IWORK( LIWORK ), LIWORK, NFO)
    call dsyevd('V','U', dimension_Hankel, hankel_tmp_eigvec, dimension_Hankel, hankel_tmp_eigval, dwork, lwork, iwork, liwork, info)
    if (info .ne. 0) then
       write(6,*)'prony_job: error! eigenvalue calculation failed for hankel_imag', info
       stop
    end if
    
    !deallocate(dwork, iwork, STAT=istat)
    
    !write(6,*)'hankel_imag_eigval'
    !write(6,1000)(dble(hankel_tmp_eigval(ni)), ni = 1,4)
    !write(6,*)
    
    !write(6,*)'hankel_imag_eigvec'
    !write(6,1000)(dble(hankel_tmp_eigvec(ni,1)), ni = 1,4)
    !write(6,1000)(dble(hankel_tmp_eigvec(ni,2)), ni = 1,4)
    !write(6,1000)(dble(hankel_tmp_eigvec(ni,3)), ni = 1,4)
    !write(6,1000)(dble(hankel_tmp_eigvec(ni,4)), ni = 1,4)
    !write(6,*)
    
    allocate(nindextmp(dimension_Hankel), STAT=istat)
    allocate(dvectmp1(dimension_Hankel), STAT=istat)
    allocate(cmattmp1(dimension_Hankel,dimension_Hankel), cmattmp2(dimension_Hankel,dimension_Hankel), STAT=istat)
    allocate(zhankel_imag_eigval(dimension_Hankel), zhankel_imag_eigvec(dimension_Hankel,dimension_Hankel), STAT=istat)
    
    nindextmp = 0
    dvectmp1 = 0.0d0
    cmattmp1 = czero
    cmattmp2 = czero
    
    do nj = 1, dimension_Hankel
        if (hankel_tmp_eigval(nj) < 0.0d0) then
            do ni = 1, dimension_Hankel
                zhankel_imag_eigvec(ni,nj) = -eye * hankel_tmp_eigvec(ni,nj)
            end do
        else if (hankel_tmp_eigval(nj) > 0.0d0) then
            do ni = 1, dimension_Hankel
                zhankel_imag_eigvec(ni,nj) = cunit * hankel_tmp_eigvec(ni,nj)
            end do
        end if
    end do
    
    do ni = 1, dimension_Hankel
        dvectmp1(ni) = dabs(hankel_tmp_eigval(ni))
        nindextmp(ni) = ni
    end do

    
    deallocate(hankel_tmp_eigval, hankel_tmp_eigvec, STAT=istat)
    
    cmattmp1 = zhankel_imag_eigvec
    
    ! sort dabs(hankel_imag_eigval) in the descending order
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
    
    zhankel_imag_eigval = dcmplx(dvectmp1, 0.0d0)
    zhankel_imag_eigvec = cmattmp2

    
    deallocate(nindextmp, STAT=istat)
    deallocate(dvectmp1, STAT=istat)
    deallocate(cmattmp1, cmattmp2, STAT=istat)
    
    !!!!!!!!!!!! check the accurancy of Takagi factorization for imaginary Hankel matirx !!!!!!!!!!!!

    allocate(cmattmp1(dimension_Hankel,dimension_Hankel), cmattmp2(dimension_Hankel,dimension_Hankel), STAT=istat)
    allocate(cmattmp3(dimension_Hankel,dimension_Hankel), cmattmp4(dimension_Hankel,dimension_Hankel), STAT=istat)
    cmattmp1 = czero
    cmattmp2 = czero
    cmattmp3 = czero
    cmattmp4 = czero
    error = 0.0d0
    
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
    call zgemm('N', 'T', dimension_Hankel, dimension_Hankel, dimension_Hankel, cunit, zhankel_imag_eigvec, &
            dimension_Hankel, cmattmp1, dimension_Hankel, czero, cmattmp2, dimension_Hankel)
          
     
    ! cmattmp3 =  zhankel_imag_eigvec * hankel_imag_eigval * zhankel_imag_eigvec.T
    call zgemm('N', 'T', dimension_Hankel, dimension_Hankel, dimension_Hankel, cunit, cmattmp2, &
            dimension_Hankel, zhankel_imag_eigvec, dimension_Hankel, czero, cmattmp3, dimension_Hankel)
    
    dtmp1 = 0.d0
    do ni = 1, dimension_Hankel
        do nj = 1, dimension_Hankel
            dtmp1 = max(dtmp1, cdabs( dcmplx(hankel_imag(ni,nj), 0.0d0) - cmattmp3(ni,nj) ))
            error = error + cdabs( dcmplx(hankel_imag(ni,nj), 0.0d0) - cmattmp3(ni,nj) )
        end do
    end do
    
    write(6,*)'prony_job: error of Takagi facotization for imagnary part', error 
    write(6,*)'prony_job: the maxium error element is', dtmp1
    write(6,*)
    
    deallocate(cmattmp1,cmattmp2,cmattmp3,cmattmp4, STAT=istat)
    
    !!!!!!!!!!!! calculate the roots of polynomial !!!!!!!!!!!!
    
    allocate(cpara(dimension_Hankel - 1), cpara_mat(dimension_Hankel - 1,dimension_Hankel - 1), STAT=istat)
    allocate(prony_roots_imag(dimension_Hankel - 1), STAT=istat)
    cpara = czero
    cpara_mat = czero
    prony_roots_imag = czero
    
    do ni = 1, dimension_Hankel - 1
        cpara(ni) = zhankel_imag_eigvec(ni, nprony_imag + 1)/(zhankel_imag_eigvec(dimension_Hankel, nprony_imag + 1))
        if (ni <= 10 .or. (dimension_Hankel - 1 - ni) <= 10) then
            write(6,*)'numerator', zhankel_imag_eigvec(ni, nprony_imag + 1)
            write(6,*)'denominator',(zhankel_imag_eigvec(dimension_Hankel, nprony_imag + 1))
        end if
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
    
    lwork = 5*(dimension_Hankel-1)
    allocate(cmattmp1(dimension_Hankel - 1, dimension_Hankel - 1),cmattmp2(dimension_Hankel - 1, dimension_Hankel - 1), STAT=istat)
    allocate(cdwork(lwork), dwork(2*(dimension_Hankel-1)), STAT=istat)
    cdwork = czero
    cmattmp1 = czero
    cmattmp2 = czero
    dwork = 0.0d0
    dtmp1 = 0.0d0
    
    !call zgeev(JOBVL, JOBVR, N, A( lda, N ), LDA, W( N ), VL( ldvl, N ), LDVL, VR( ldvr, N ), LDVR, WORK( LWORK ), LWORK, RWORK( 2N ), INFO)
    call zgeev('N','N', dimension_Hankel - 1, cpara_mat, dimension_Hankel - 1, prony_roots_imag, &
            cmattmp1, dimension_Hankel - 1, cmattmp2, dimension_Hankel - 1, cdwork, 5*(dimension_Hankel-1), dwork, info)
    
    deallocate(cdwork, STAT=istat)
    deallocate(cmattmp1,cmattmp2, STAT=istat)
    deallocate(zhankel_imag_eigvec, zhankel_imag_eigval, STAT=istat)
    
    !!!!!!!!!!!! get the number of nonzero roots for imagniary part !!!!!!!!!!!!
    
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
    
        write(6,*) 'prony_job: there is no zero element'
        write(6,*) 'prony_job: number of nonzero element', ncount
        write(6,*)
    else if (ncount .ne. (dimension_Hankel - 1)) then
        deallocate(prony_roots_imag)
        allocate(prony_roots_imag(ncount), STAT=istat)
        write(6,*) 'prony_job: number of nonzero element', ncount
        write(6,*)
        do ni = 1, ncount
            prony_roots_imag(ni) = cvectmp1(ni)
        end do
    end if
    
    deallocate(cvectmp1, STAT=istat)
        
    !!!!!!!!!!!! get the order of dabs(prony_roots_imag) in dabs(prony_roots_imag) sorted in the ascending order !!!!!!!!!!!! 
       
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
                
                ntmp1 = nindextmp(nj)
                nindextmp(nj) = nindextmp(nj + 1)
                nindextmp(nj + 1) = ntmp1
            end if
         end do
    end do
    
    allocate(cvectmp1(nprony_imag), STAT=istat)
    cvectmp1 = czero
    
    do ni = 1, nprony_imag
        cvectmp1(ni) = prony_roots_imag(nindextmp(ni))
    end do
    
    deallocate(prony_roots_imag, STAT=istat)
    
    !!!!!!!!!!!! calculate the dissipation rates of imagniary part !!!!!!!!!!!!
    
    allocate(prony_roots_imag(nprony_imag),diss_rate_prony_imag(nprony_imag), STAT=istat)
    prony_roots_imag = cvectmp1
    diss_rate_prony_imag = czero
    
    cvectmp1 = czero
    do ni = 1, nprony_imag
        ctmp1 = dcmplx(atan( dimag(prony_roots_imag(ni)), dble(prony_roots_imag(ni))), 0.0d0)
        cvectmp1(ni) = -2.d0 * (dimension_Hankel - 1) / tt_for_prony(tt_sample_for_prony) &
                  * ( dlog( cdabs(prony_roots_imag(ni)) ) + eye * ctmp1)
    end do
    
    diss_rate_prony_imag = cvectmp1
    deallocate(cvectmp1, STAT=istat)
    
    !write(6,*)'diss_rate_prony_imag'
    !write(6,*) (diss_rate_prony_imag(ni), ni = 1, nprony)
    !write(6,*)
    
    !!!!!!!!!!!! calculate the coefficients of imagniary part !!!!!!!!!!!!
    ! calculate the coefficients by the least squares method
    
    allocate(cmattmp1(tt_sample_for_prony, nprony_imag), STAT=istat)
    allocate(cmattmp2(nprony_imag, tt_sample_for_prony), STAT=istat)
    allocate(cmattmp3(nprony_imag, nprony_imag), STAT=istat)
    allocate(cmattmp4(nprony_imag, nprony_imag), STAT=istat)
    allocate(cvectmp1(tt_sample_for_prony), STAT=istat)
    allocate(nvectmp1(nprony_imag), STAT=istat)
    cmattmp1 = czero
    cmattmp2 = czero
    cmattmp3 = czero
    cmattmp4 = czero
    cvectmp1 = czero
    nvectmp1 = 0
    
    do ni = 1, tt_sample_for_prony
        do nj = 1, nprony_imag
            cmattmp1(ni,nj) = prony_roots_imag(nj)**(ni - 1)
            cmattmp2(nj,ni) = prony_roots_imag(nj)**(ni - 1)
        end do
    end do
    
    ! cmattmp3 = cmattmp2 * cmattmp1 = omega_imag.T * omega_imag
    call zgemm('N', 'N', nprony_imag, nprony_imag, tt_sample_for_prony, cunit, cmattmp2, &
             nprony_imag, cmattmp1, tt_sample_for_prony, czero, cmattmp3, nprony_imag)
    
    ! cvectmp1 = cmattmp2 * dimag(res_prony) = omega_imag.T * dimag(res_prony)
    do ni = 1, nprony_imag
        do nj = 1, tt_sample_for_prony
            cvectmp1(ni) = cvectmp1(ni) + cmattmp2(ni,nj) * dimag(res_prony(nj))
        end do
    end do
    
    ! inverse cmattmp3
    lwork = 5 * nprony_imag
    allocate(cvectmp2(lwork), STAT=istat)
    cvectmp2 = czero
    cmattmp4 = cmattmp3
    
    ! calculate the coefficients by the least squares method
    call zgetrf(nprony_imag, nprony_imag, cmattmp3, nprony_imag, nvectmp1, info)
    call zgetri(nprony_imag, cmattmp3, nprony_imag, nvectmp1, cvectmp2, lwork, info) 
    
    deallocate(nvectmp1, cvectmp2, STAT=istat)
    
    
    allocate(cmattmp5(nprony_imag, nprony_imag), STAT=istat)
    
    cmattmp5 = czero
    call zgemm('N', 'N', nprony_imag, nprony_imag, nprony_imag, cunit, cmattmp4, &
             nprony_imag, cmattmp3, nprony_imag, czero, cmattmp5, nprony_imag)
    
    allocate(eta_prony_imag(nprony_imag), STAT=istat)
    eta_prony_imag = czero
    
    ! cmattmp3 * cvectmp1
    do ni = 1, nprony_imag
        do nj = 1, nprony_imag
            eta_prony_imag(ni) = eta_prony_imag(ni) + cmattmp3(ni,nj) * cvectmp1(nj)
        end do
    end do
    
    deallocate(cmattmp1, cmattmp2, cmattmp3, cmattmp4, cmattmp5, cvectmp1, STAT=istat)
    
    !write(6,*)'eta_prony_imag'
    !write(6,*) (eta_prony_imag(ni), ni = 1, nprony)
    !write(6,*)
    
    
    !!!!!!!!!!!! calculate the dissipation rates of real part !!!!!!!!!!!!
    
    allocate(diss_rate_prony_real(nprony_real), STAT=istat)
    diss_rate_prony_real(1) = band_width / hbar

    allocate(eta_prony_real(nprony_real), STAT=istat)
    eta_prony_real(1) = (sys_bath_coupling * band_width * 0.25d0) / (hbar**2)
    
    !write(6,*)'eta_prony_real'
    !write(6,*) (eta_prony_real(ni), ni = 1, nprony)
    !write(6,*)
    
    !!!!!!!!!!!! output dissipation rate and coefficients !!!!!!!!!!!!
    
    allocate(eta_prony(nprony_imag + nprony_imag), diss_rate_prony(nprony_imag + nprony_imag), STAT=istat)
    do ni = 1, (nprony_imag + nprony_real)
        if (ni <= nprony_imag) then
            eta_prony(ni) = eta_prony_imag(ni) * eye
            diss_rate_prony(ni) = diss_rate_prony_imag(ni)
        else if (ni > nprony_imag) then
            eta_prony(ni) = eta_prony_real(ni - nprony_imag)
            diss_rate_prony(ni) = diss_rate_prony_real(ni - nprony_imag)
        end if
    end do
    
    do ni = 1, (nprony_imag + nprony_real)
        dtmp1 = dble(eta_prony(ni))
        dtmp2 = dimag(eta_prony(ni))
        if (dabs(dtmp1) <= 1.d-6) then
            dtmp1 = 0.d0
        end if
        if (dabs(dtmp2) <= 1.d-6) then
            dtmp2 = 0.d0
        end if
        eta_prony(ni) = dcmplx(dtmp1,dtmp2)
        
        dtmp1 = dble(diss_rate_prony(ni))
        dtmp2 = dimag(diss_rate_prony(ni))
        if (dabs(dtmp1) <= 1.d-6) then
            dtmp1 = 0.d0
        end if
        if (dabs(dtmp2) <= 1.d-6) then
            dtmp2 = 0.d0
        end if
        diss_rate_prony(ni) = dcmplx(dtmp1,dtmp2)
    end do
    
    open(unit = 90, file='eta_real.data')
    rewind(90)
    open(unit = 91, file='eta_imag.data')
    rewind(91)
    open(unit = 92, file='rate_real.data')
    rewind(92)
    open(unit = 93, file='rate_imag.data')
    rewind(93)
    
    do ni =1, (nprony_imag + nprony_real)
        write(90,*) dble(eta_prony(ni))
        write(91,*) dimag(eta_prony(ni))
        write(92,*) dble(diss_rate_prony(ni))
        write(93,*) dimag(diss_rate_prony(ni))
    end do
    
    close(90)
    close(91)
    close(92)
    close(93)
    
    ! only valid for bcenter = 0
    open(unit = 94, file='cb.data')
    rewind(94)
    open(unit = 95, file='cd.data')
    rewind(95)
    open(unit = 96, file='cgama.data')
    rewind(96)
    write(94,*)'isgn, ialf, icor, ispin, iorbs2, Re[cb], Imag[cb]'
    write(95,*)'isgn, ialf, icor, ispin, iorbs2, Re[cd], Imag[cd]'
    write(96,*)'isgn, ialf, icor, ispin, iorbs2, Re[cgama], Imag[cgama]'
    do nm = 1, 2
        do nl = 1, 1
            do nk = 1, (nprony_real + nprony_imag)
                do nj = 1, 2
                    do ni = 1, 1
                        !if (nm == 1) then
                            write(94,1001)nm, nl, nk, nj, ni, dble(eta_prony(nk)), dimag(eta_prony(nk))
                            write(95,1001)nm, nl, nk, nj, ni, dble(eta_prony(nk)), -dimag(eta_prony(nk))
                            write(96,1001)nm, nl, nk, nj, ni, -dble(diss_rate_prony(nk)), -dimag(diss_rate_prony(nk))
                        !else if (nm == 2) then
                        !    write(94,1001)nm, nl, nk, nj, ni, dble(eta_prony(nk)), dimag(eta_prony(nk))
                        !    write(95,1001)nm, nl, nk, nj, ni, dble(eta_prony(nk)), -dimag(eta_prony(nk))
                        !    write(96,1001)nm, nl, nk, nj, ni, dble(diss_rate_prony(nk)), -dimag(diss_rate_prony(nk))
                        !end if
                    end do
                end do
            end do
        end do
    end do
    close(94)
    close(95)
    close(96)
    
    ! only valid for bcenter = 0
    !open(unit = 94, file='cb.data')
    !rewind(94)
    !open(unit = 95, file='cd.data')
    !rewind(95)
    !open(unit = 96, file='cgama.data')
    !rewind(96)
    !write(94,*)'Re[cb], Imag[cb]'
    !write(95,*)'Re[cd], Imag[cd]'
    !write(96,*)'Re[cgama], Imag[cgama]'
    !do nm = 1, 2
    !    do nl = 1, 1
    !        do nk = 1, (nprony_real + nprony_imag)
    !            do nj = 1, 2
    !                do ni = 1, 1
    !                    write(94,*)dble(eta_prony(nk)), dimag(eta_prony(nk))
    !                    write(95,*)dble(eta_prony(nk)), -dimag(eta_prony(nk))
    !                    write(96,*)-dble(diss_rate_prony(nk)), -dimag(diss_rate_prony(nk))
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !close(94)
    !close(95)
    !close(96)

     
    write(6,*)'eta_prony'
    write(6,*) (eta_prony(ni), ni = 1, (nprony_real + nprony_imag))    
    write(6,*)'diss_rate_prony'
    write(6,*) (diss_rate_prony(ni), ni = 1, (nprony_real + nprony_imag))
    write(6,*)
    
    !do ni =1, (nprony_real + nprony_imag)
    !    output_eta(ni) = eta_prony(ni)
    !    output_diss_rate(ni) = diss_rate_prony(ni)
    !end do
    
    1001 format(5(I3,1x),2(1x, f37.32))
    1000 format(20(f17.13,2x))
end program prony_job
!end subroutine prony_job
