program read_diff_check
implicit none
include "../../include/sizes"
include "../../include/common"
integer :: nalf_test
!
! 独立测试程序：用于单独检查 input 中 $diffdecomp 和后续 scheme 行。
! 放在 tools/read_diff_check 下，避免 main/*.f90 自动编译时产生第二个主程序。
!
nalf_test = 2
open(unit=5,file='input',status='old')
call read_diffdecomp_config(5, nalf_test)
close(5)
write(6,*)'read_diff_check: diffdecomp parser finished.'
end program read_diff_check
