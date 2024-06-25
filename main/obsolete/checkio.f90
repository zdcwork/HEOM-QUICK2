subroutine checkio
implicit none
include '../include/sizes'
include '../include/common'
!
write(6,*)
write(6,*)'checkio: check the status of all asynchronous files '
call flush(6)
!
write(6,*)'checkio: file unit=17 '
call flush(6)
wait(unit=17)
close(17)
write(6,*)'checkio: file unit=17 closed '
call flush(6)
!
return
!
end subroutine checkio
