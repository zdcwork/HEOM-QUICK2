subroutine zpout_psf(mpole, msgn, zcoef, zpole)
!
! purpose : find poles and coefficients for PSF scheme
!
implicit none
include "../include/sizes"
include "../include/common"
!
integer, intent (in)     :: mpole, msgn
complex*16, intent (out) :: zcoef(*), zpole(mpole,*)
integer                  :: isgn
!
if (.not. psfjob) then
    write(6,*)'zpout_psf: error! psfjob= ', psfjob
    stop
end if
!
if (itype_psf .eq. 1) then
    if (mpole .ne. npsf1) then
        write(6,*)'zpout_psf: error! mpole != npsf1 found '
        write(6,*)mpole, npsf1
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf1(1:mpole), coefi_psf1(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf1(1:mpole), gamai_psf1(1:mpole))
    end do
else if (itype_psf .eq. 2) then
    if (mpole .ne. npsf2) then
        write(6,*)'zpout_psf: error! mpole != npsf2 found '
        write(6,*)mpole, npsf2
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf2(1:mpole), coefi_psf2(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf2(1:mpole), gamai_psf2(1:mpole))
    end do
else if (itype_psf .eq. 3) then
    if (mpole .ne. npsf3) then
        write(6,*)'zpout_psf: error! mpole != npsf3 found '
        write(6,*)mpole, npsf3
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf3(1:mpole), coefi_psf3(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf3(1:mpole), gamai_psf3(1:mpole))
    end do
else if (itype_psf .eq. 4) then
    if (mpole .ne. npsf4) then
        write(6,*)'zpout_psf: error! mpole != npsf4 found '
        write(6,*)mpole, npsf4
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf4(1:mpole), coefi_psf4(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf4(1:mpole), gamai_psf4(1:mpole))
    end do
else if (itype_psf .eq. 5) then
    if (mpole .ne. npsf5) then
        write(6,*)'zpout_psf: error! mpole != npsf5 found '
        write(6,*)mpole, npsf5
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf5(1:mpole), coefi_psf5(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf5(1:mpole), gamai_psf5(1:mpole))
    end do
else if (itype_psf .eq. 6) then
    if (mpole .ne. npsf6) then
        write(6,*)'zpout_psf: error! mpole != npsf6 found '
        write(6,*)mpole, npsf6
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf6(1:mpole), coefi_psf6(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf6(1:mpole), gamai_psf6(1:mpole))
    end do
else if (itype_psf .eq. 7) then
    if (mpole .ne. npsf7) then
        write(6,*)'zpout_psf: error! mpole != npsf7 found '
        write(6,*)mpole, npsf7
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf7(1:mpole), coefi_psf7(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf7(1:mpole), gamai_psf7(1:mpole))
    end do
else if (itype_psf .eq. 8) then
    if (mpole .ne. npsf8) then
        write(6,*)'zpout_psf: error! mpole != npsf8 found '
        write(6,*)mpole, npsf8
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf8(1:mpole), coefi_psf8(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf8(1:mpole), gamai_psf8(1:mpole))
    end do
else if (itype_psf .eq. 9) then
    if (mpole .ne. npsf9) then
        write(6,*)'zpout_psf: error! mpole != npsf9 found '
        write(6,*)mpole, npsf9
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf9(1:mpole), coefi_psf9(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf9(1:mpole), gamai_psf9(1:mpole))
    end do
else if (itype_psf .eq. 10) then
    if (mpole .ne. npsf10) then
        write(6,*)'zpout_psf: error! mpole != npsf10 found '
        write(6,*)mpole, npsf10
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf10(1:mpole), coefi_psf10(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf10(1:mpole), gamai_psf10(1:mpole))
    end do
else if (itype_psf .eq. 11) then
    if (mpole .ne. npsf11) then
        write(6,*)'zpout_psf: error! mpole != npsf11 found '
        write(6,*)mpole, npsf11
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf11(1:mpole), coefi_psf11(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf11(1:mpole), gamai_psf11(1:mpole))
    end do
else if (itype_psf .eq. 12) then
    if (mpole .ne. npsf12) then
        write(6,*)'zpout_psf: error! mpole != npsf12 found '
        write(6,*)mpole, npsf12
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf12(1:mpole), coefi_psf12(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf12(1:mpole), gamai_psf12(1:mpole))
    end do
else if (itype_psf .eq. 13) then
    if (mpole .ne. npsf13) then
        write(6,*)'zpout_psf: error! mpole != npsf13 found '
        write(6,*)mpole, npsf13
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf13(1:mpole), coefi_psf13(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf13(1:mpole), gamai_psf13(1:mpole))
    end do
else if (itype_psf .eq. 14) then
    if (mpole .ne. npsf14) then
        write(6,*)'zpout_psf: error! mpole != npsf14 found '
        write(6,*)mpole, npsf14
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf14(1:mpole), coefi_psf14(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf14(1:mpole), gamai_psf14(1:mpole))
    end do
else if (itype_psf .eq. 15) then
    if (mpole .ne. npsf15) then
        write(6,*)'zpout_psf: error! mpole != npsf15 found '
        write(6,*)mpole, npsf15
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf15(1:mpole), coefi_psf15(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf15(1:mpole), gamai_psf15(1:mpole))
    end do
else if (itype_psf .eq. 16) then
    if (mpole .ne. npsf16) then
        write(6,*)'zpout_psf: error! mpole != npsf16 found '
        write(6,*)mpole, npsf16
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf16(1:mpole), coefi_psf16(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf16(1:mpole), gamai_psf16(1:mpole))
    end do
else if (itype_psf .eq. 17) then
    if (mpole .ne. npsf17) then
        write(6,*)'zpout_psf: error! mpole != npsf17 found '
        write(6,*)mpole, npsf17
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf17(1:mpole), coefi_psf17(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf17(1:mpole), gamai_psf17(1:mpole))
    end do
else if (itype_psf .eq. 18) then
    if (mpole .ne. npsf18) then
        write(6,*)'zpout_psf: error! mpole != npsf18 found '
        write(6,*)mpole, npsf18
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf18(1:mpole), coefi_psf18(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf18(1:mpole), gamai_psf18(1:mpole))
    end do
else if (itype_psf .eq. 19) then
    if (mpole .ne. npsf19) then
        write(6,*)'zpout_psf: error! mpole != npsf19 found '
        write(6,*)mpole, npsf19
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf19(1:mpole), coefi_psf19(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf19(1:mpole), gamai_psf19(1:mpole))
    end do
else if (itype_psf .eq. 20) then
    if (mpole .ne. npsf20) then
        write(6,*)'zpout_psf: error! mpole != npsf20 found '
        write(6,*)mpole, npsf20
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf20(1:mpole), coefi_psf20(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf20(1:mpole), gamai_psf20(1:mpole))
    end do
else if (itype_psf .eq. 21) then
    if (mpole .ne. npsf21) then
        write(6,*)'zpout_psf: error! mpole != npsf21 found '
        write(6,*)mpole, npsf21
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf21(1:mpole), coefi_psf21(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf21(1:mpole), gamai_psf21(1:mpole))
    end do
else if (itype_psf .eq. 22) then
    if (mpole .ne. npsf22) then
        write(6,*)'zpout_psf: error! mpole != npsf22 found '
        write(6,*)mpole, npsf22
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf22(1:mpole), coefi_psf22(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf22(1:mpole), gamai_psf22(1:mpole))
    end do
else if (itype_psf .eq. 23) then
    if (mpole .ne. npsf23) then
        write(6,*)'zpout_psf: error! mpole != npsf23 found '
        write(6,*)mpole, npsf23
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf23(1:mpole), coefi_psf23(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf23(1:mpole), gamai_psf23(1:mpole))
    end do
else if (itype_psf .eq. 24) then
    if (mpole .ne. npsf24) then
        write(6,*)'zpout_psf: error! mpole != npsf24 found '
        write(6,*)mpole, npsf24
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf24(1:mpole), coefi_psf24(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf24(1:mpole), gamai_psf24(1:mpole))
    end do
else if (itype_psf .eq. 25) then
    if (mpole .ne. npsf25) then
        write(6,*)'zpout_psf: error! mpole != npsf25 found '
        write(6,*)mpole, npsf25
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf25(1:mpole), coefi_psf25(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf25(1:mpole), gamai_psf25(1:mpole))
    end do
else if (itype_psf .eq. 26) then
    if (mpole .ne. npsf26) then
        write(6,*)'zpout_psf: error! mpole != npsf26 found '
        write(6,*)mpole, npsf26
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf26(1:mpole), coefi_psf26(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf26(1:mpole), gamai_psf26(1:mpole))
    end do
else if (itype_psf .eq. 27) then
    if (mpole .ne. npsf27) then
        write(6,*)'zpout_psf: error! mpole != npsf27 found '
        write(6,*)mpole, npsf27
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf27(1:mpole), coefi_psf27(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf27(1:mpole), gamai_psf27(1:mpole))
    end do
else if (itype_psf .eq. 28) then
    if (mpole .ne. npsf28) then
        write(6,*)'zpout_psf: error! mpole != npsf28 found '
        write(6,*)mpole, npsf28
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf28(1:mpole), coefi_psf28(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf28(1:mpole), gamai_psf28(1:mpole))
    end do
else if (itype_psf .eq. 29) then
    if (mpole .ne. npsf29) then
        write(6,*)'zpout_psf: error! mpole != npsf29 found '
        write(6,*)mpole, npsf29
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf29(1:mpole), coefi_psf29(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf29(1:mpole), gamai_psf29(1:mpole))
    end do
else if (itype_psf .eq. 30) then
    if (mpole .ne. npsf30) then
        write(6,*)'zpout_psf: error! mpole != npsf30 found '
        write(6,*)mpole, npsf30
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf30(1:mpole), coefi_psf30(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf30(1:mpole), gamai_psf30(1:mpole))
    end do
else if (itype_psf .eq. 31) then
    if (mpole .ne. npsf31) then
        write(6,*)'zpout_psf: error! mpole != npsf31 found '
        write(6,*)mpole, npsf31
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf31(1:mpole), coefi_psf31(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf31(1:mpole), gamai_psf31(1:mpole))
    end do
else if (itype_psf .eq. 32) then
    if (mpole .ne. npsf32) then
        write(6,*)'zpout_psf: error! mpole != npsf32 found '
        write(6,*)mpole, npsf32
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf32(1:mpole), coefi_psf32(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf32(1:mpole), gamai_psf32(1:mpole))
    end do
else if (itype_psf .eq. 33) then
    if (mpole .ne. npsf33) then
        write(6,*)'zpout_psf: error! mpole != npsf33 found '
        write(6,*)mpole, npsf33
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf33(1:mpole), coefi_psf33(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf33(1:mpole), gamai_psf33(1:mpole))
    end do
else if (itype_psf .eq. 34) then
    if (mpole .ne. npsf34) then
        write(6,*)'zpout_psf: error! mpole != npsf34 found '
        write(6,*)mpole, npsf34
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf34(1:mpole), coefi_psf34(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf34(1:mpole), gamai_psf34(1:mpole))
    end do
else if (itype_psf .eq. 35) then
    if (mpole .ne. npsf35) then
        write(6,*)'zpout_psf: error! mpole != npsf35 found '
        write(6,*)mpole, npsf35
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf35(1:mpole), coefi_psf35(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf35(1:mpole), gamai_psf35(1:mpole))
    end do
else if (itype_psf .eq. 36) then
    if (mpole .ne. npsf36) then
        write(6,*)'zpout_psf: error! mpole != npsf36 found '
        write(6,*)mpole, npsf36
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf36(1:mpole), coefi_psf36(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf36(1:mpole), gamai_psf36(1:mpole))
    end do
else if (itype_psf .eq. 37) then
    if (mpole .ne. npsf37) then
        write(6,*)'zpout_psf: error! mpole != npsf37 found '
        write(6,*)mpole, npsf37
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf37(1:mpole), coefi_psf37(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf37(1:mpole), gamai_psf37(1:mpole))
    end do
else if (itype_psf .eq. 38) then
    if (mpole .ne. npsf38) then
        write(6,*)'zpout_psf: error! mpole != npsf38 found '
        write(6,*)mpole, npsf38
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf38(1:mpole), coefi_psf38(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf38(1:mpole), gamai_psf38(1:mpole))
    end do
else if (itype_psf .eq. 39) then
    if (mpole .ne. npsf39) then
        write(6,*)'zpout_psf: error! mpole != npsf39 found '
        write(6,*)mpole, npsf39
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf39(1:mpole), coefi_psf39(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf39(1:mpole), gamai_psf39(1:mpole))
    end do
else if (itype_psf .eq. 40) then
    if (mpole .ne. npsf40) then
        write(6,*)'zpout_psf: error! mpole != npsf40 found '
        write(6,*)mpole, npsf40
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf40(1:mpole), coefi_psf40(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf40(1:mpole), gamai_psf40(1:mpole))
    end do
else if (itype_psf .eq. 41) then
    if (mpole .ne. npsf41) then
        write(6,*)'zpout_psf: error! mpole != npsf41 found '
        write(6,*)mpole, npsf41
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf41(1:mpole), coefi_psf41(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf41(1:mpole), gamai_psf41(1:mpole))
    end do
else if (itype_psf .eq. 42) then
    if (mpole .ne. npsf42) then
        write(6,*)'zpout_psf: error! mpole != npsf42 found '
        write(6,*)mpole, npsf42
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf42(1:mpole), coefi_psf42(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf42(1:mpole), gamai_psf42(1:mpole))
    end do
else if (itype_psf .eq. 43) then
    if (mpole .ne. npsf43) then
        write(6,*)'zpout_psf: error! mpole != npsf43 found '
        write(6,*)mpole, npsf43
        stop
    end if
    zcoef(1:mpole) = dcmplx(coefr_psf43(1:mpole), coefi_psf43(1:mpole))
    do isgn=1,msgn
       zpole(1:mpole,isgn) = dcmplx(gamar_psf43(1:mpole), gamai_psf43(1:mpole))
    end do
else
    write(6,*)'zpout_psf: error! unknown itype_psf= ', itype_psf
end if
!
return
!
end subroutine zpout_psf
