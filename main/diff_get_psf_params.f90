subroutine diff_get_psf_params(itype, npole, rpsf, tpsf)
implicit none
include "../include/sizes"
!
! 根据 PSF 类型返回预设 pole 数和参考参数。
! diff 模式中每个电极可以有自己的 itype_psf，因此不能再使用原程序
! evaluatepara.f90 中会覆盖全局 nmats/bandwidth/tpsf 的长分支。
!
integer, intent(in)  :: itype
integer, intent(out) :: npole
real*8,  intent(out) :: rpsf, tpsf
!
select case (itype)
case (1);  npole = npsf1;  rpsf = rpsf1;  tpsf = tpsf1
case (2);  npole = npsf2;  rpsf = rpsf2;  tpsf = tpsf2
case (3);  npole = npsf3;  rpsf = rpsf3;  tpsf = tpsf3
case (4);  npole = npsf4;  rpsf = rpsf4;  tpsf = tpsf4
case (5);  npole = npsf5;  rpsf = rpsf5;  tpsf = tpsf5
case (6);  npole = npsf6;  rpsf = rpsf6;  tpsf = tpsf6
case (7);  npole = npsf7;  rpsf = rpsf7;  tpsf = tpsf7
case (8);  npole = npsf8;  rpsf = rpsf8;  tpsf = tpsf8
case (9);  npole = npsf9;  rpsf = rpsf9;  tpsf = tpsf9
case (10); npole = npsf10; rpsf = rpsf10; tpsf = tpsf10
case (11); npole = npsf11; rpsf = rpsf11; tpsf = tpsf11
case (12); npole = npsf12; rpsf = rpsf12; tpsf = tpsf12
case (13); npole = npsf13; rpsf = rpsf13; tpsf = tpsf13
case (14); npole = npsf14; rpsf = rpsf14; tpsf = tpsf14
case (15); npole = npsf15; rpsf = rpsf15; tpsf = tpsf15
case (16); npole = npsf16; rpsf = rpsf16; tpsf = tpsf16
case (17); npole = npsf17; rpsf = rpsf17; tpsf = tpsf17
case (18); npole = npsf18; rpsf = rpsf18; tpsf = tpsf18
case (19); npole = npsf19; rpsf = rpsf19; tpsf = tpsf19
case (20); npole = npsf20; rpsf = rpsf20; tpsf = tpsf20
case (21); npole = npsf21; rpsf = rpsf21; tpsf = tpsf21
case (22); npole = npsf22; rpsf = rpsf22; tpsf = tpsf22
case (23); npole = npsf23; rpsf = rpsf23; tpsf = tpsf23
case (24); npole = npsf24; rpsf = rpsf24; tpsf = tpsf24
case (25); npole = npsf25; rpsf = rpsf25; tpsf = tpsf25
case (26); npole = npsf26; rpsf = rpsf26; tpsf = tpsf26
case (27); npole = npsf27; rpsf = rpsf27; tpsf = tpsf27
case (28); npole = npsf28; rpsf = rpsf28; tpsf = tpsf28
case (29); npole = npsf29; rpsf = rpsf29; tpsf = tpsf29
case (30); npole = npsf30; rpsf = rpsf30; tpsf = tpsf30
case (31); npole = npsf31; rpsf = rpsf31; tpsf = tpsf31
case (32); npole = npsf32; rpsf = rpsf32; tpsf = tpsf32
case (33); npole = npsf33; rpsf = rpsf33; tpsf = tpsf33
case (34); npole = npsf34; rpsf = rpsf34; tpsf = tpsf34
case (35); npole = npsf35; rpsf = rpsf35; tpsf = tpsf35
case (36); npole = npsf36; rpsf = rpsf36; tpsf = tpsf36
case (37); npole = npsf37; rpsf = rpsf37; tpsf = tpsf37
case (38); npole = npsf38; rpsf = rpsf38; tpsf = tpsf38
case (39); npole = npsf39; rpsf = rpsf39; tpsf = tpsf39
case (40); npole = npsf40; rpsf = rpsf40; tpsf = tpsf40
case (41); npole = npsf41; rpsf = rpsf41; tpsf = tpsf41
case (42); npole = npsf42; rpsf = rpsf42; tpsf = tpsf42
case (43); npole = npsf43; rpsf = rpsf43; tpsf = tpsf43
case default
   write(6,*)'diff_get_psf_params: error! unsupported itype_psf = ', itype
   stop
end select
end subroutine diff_get_psf_params
