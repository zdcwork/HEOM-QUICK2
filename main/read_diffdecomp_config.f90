subroutine read_diffdecomp_config(iunit, nalf_in)
implicit none
include "../include/sizes"
include "../include/common"
!
! 读取不同电极谱分解配置。
! iunit    : 已打开的 input 文件单元，主程序中通常为 5。
! nalf_in  : 当前体系的电极数；配置行必须覆盖 1:nalf_in 的每个电极。
!
integer, intent(in) :: iunit, nalf_in
logical             :: seen(maxalf)
integer             :: ialf, ischeme, ipoles, itype
integer             :: ios, iline, ndict
character(len=256)  :: line
character(len=64)   :: scheme
character(len=64)   :: dictionary(9)
namelist / diffdecomp / ldiff_spe_dec, ldiff_spe_dec_sp
!
dictionary(1) = 'pade'
dictionary(2) = 'psdjob'
dictionary(3) = 'prony'
dictionary(4) = 'psfjob'
dictionary(5) = 'psdfff'
dictionary(6) = 'fano'
dictionary(7) = 'read'
dictionary(8) = 'mtsjob'
dictionary(9) = 'lprony'
ndict = 9
!
! 默认关闭 diff 模式，并清空所有按电极配置。
ldiff_spe_dec    = .false.
ldiff_spe_dec_sp = .false.
ncount_alf       = 0
ncount_scheme    = 0
ncount_poles     = 0
ncount_type      = 0
nmats_alf        = 0
ncor_alf         = 0
lpsdfff_alf      = .false.
seen             = .false.
!
rewind(iunit)
read(iunit, nml=diffdecomp, iostat=ios, end=100)
if (ios .ne. 0) then
   write(6,*)'read_diffdecomp_config: error reading $diffdecomp namelist, iostat = ', ios
   stop
end if
100 continue
!
! 未启用 diff 模式时立即返回，保持原程序行为。
if (.not. (ldiff_spe_dec .or. ldiff_spe_dec_sp)) return
!
write(6,*)
write(6,*)'Different spectral decomposition schemes are used for different leads.'
!
do iline = 1, nalf_in
   read(iunit,'(A)',iostat=ios) line
   if (ios .ne. 0) then
      write(6,*)'read_diffdecomp_config: error! not enough scheme lines after $diffdecomp'
      write(6,*)'read_diffdecomp_config: expected nalf lines, nalf = ', nalf_in
      stop
   end if
!
   ialf   = 0
   scheme = ' '
   ipoles = 0
   itype  = 0
!
! 先尝试读取可选 type；若该列不存在，再按三列格式读取。
   read(line,*,iostat=ios) ialf, scheme, ipoles, itype
   if (ios .ne. 0) then
      itype = 0
      read(line,*,iostat=ios) ialf, scheme, ipoles
   end if
   if (ios .ne. 0) then
      write(6,*)'read_diffdecomp_config: error parsing scheme line: ', trim(line)
      stop
   end if
!
! 检查电极编号范围和重复指定。
   if (ialf .lt. 1 .or. ialf .gt. nalf_in) then
      write(6,*)'read_diffdecomp_config: error! ialf out of range: ', ialf
      write(6,*)'read_diffdecomp_config: valid range is 1:nalf = ', nalf_in
      stop
   end if
   if (seen(ialf)) then
      write(6,*)'read_diffdecomp_config: error! duplicated scheme for ialf = ', ialf
      stop
   end if
   seen(ialf) = .true.
!
! scheme 名称大小写不敏感，统一转成小写再查字典。
   call diff_to_lower(scheme)
   ischeme = 0
   do ios = 1, ndict
      if (trim(scheme) .eq. trim(dictionary(ios))) then
         ischeme = ios
         exit
      end if
   end do
   if (ischeme .eq. 0) then
      write(6,*)'read_diffdecomp_config: error! unknown scheme name: ', trim(scheme)
      stop
   end if
!
! pole 数必须为正。read/lprony 分支也使用该值作为期望/最小配置项数。
   if (ipoles .le. 0) then
      write(6,*)'read_diffdecomp_config: error! pole number must be positive for ialf = ', ialf
      write(6,*)'read_diffdecomp_config: pole number = ', ipoles
      stop
   end if
!
! Fano 和 PSF 需要 type 参数；此处只检查是否给出正整数，具体范围在 evaluatepara_diff 中进一步检查。
   if ((ischeme .eq. 3 .or. ischeme .eq. 4 .or. ischeme .eq. 5 .or. ischeme .eq. 6) .and. itype .le. 0) then
      write(6,*)'read_diffdecomp_config: error! scheme needs a positive type parameter.'
      write(6,*)'read_diffdecomp_config: ialf, scheme, type = ', ialf, trim(scheme), itype
      stop
   end if
!
   ncount_alf(ialf)    = ialf
   ncount_scheme(ialf) = ischeme
   ncount_poles(ialf)  = ipoles
   ncount_type(ialf)   = itype
!
   if (ischeme .eq. 1 .or. ischeme .eq. 2) then
      write(6,1000) ialf, 'Pade', ipoles
   else if (ischeme .eq. 3 .or. ischeme .eq. 4) then
      write(6,1001) ialf, 'PSF', ipoles, itype
   else if (ischeme .eq. 5 .or. ischeme .eq. 6) then
      write(6,1001) ialf, 'Fano', ipoles, itype
   else if (ischeme .eq. 8) then
      write(6,1000) ialf, 'Matsubara', ipoles
   else
      write(6,1000) ialf, 'external Prony/read', ipoles
   end if
end do
!
do ialf = 1, nalf_in
   if (.not. seen(ialf)) then
      write(6,*)'read_diffdecomp_config: error! no scheme specified for ialf = ', ialf
      stop
   end if
end do
call flush(6)
return
!
1000 format('For lead ', I4, ' using ', A, ' scheme with nmats = ', I6)
1001 format('For lead ', I4, ' using ', A, ' scheme with nmats = ', I6, ' and itype = ', I6)
end subroutine read_diffdecomp_config
!
subroutine diff_to_lower(str)
implicit none
character(len=*), intent(inout) :: str
integer :: i, ia, iz, ishift, ic
!
ia = ichar('A')
iz = ichar('Z')
ishift = ichar('a') - ichar('A')
do i = 1, len_trim(str)
   ic = ichar(str(i:i))
   if (ic .ge. ia .and. ic .le. iz) str(i:i) = char(ic + ishift)
end do
end subroutine diff_to_lower
