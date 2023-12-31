subroutine genoutput
implicit none
include '../include/sizes'
include '../include/common'
!
character*10            :: date, time
character*120           :: line
integer                 :: istat
!
call date_and_time(date, time)
!
write(6,*)
write(6,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
write(6,*)
write(6,*)'                      H E O M - Q U I C K                                '
write(6,*)
Write(6,*)'              Hierarchical Equations of Motion for                       '
write(6,*)'            QUantum Impurity with a Correlated Kernel                    '
write(6,*)
write(6,*)'                          developed by                                   '
write(6,*)'         Prof. Xiao Zheng and Prof. YiJing Yan''s Groups @               ' 
write(6,*)'          University of Science and Technology of China                  '
write(6,*)'                        Fudan University                                 '    
write(6,*)
write(6,*)'            Email: xz58@ustc.edu.cn, xzheng@fudan.edu.cn                 '
write(6,*)
write(6,1001)version
write(6,*)
write(6,*)'           CopyRight (C) 2007-2023. All Rights Reserved.                 '
write(6,*)
write(6,*)'  References:                                                            '
write(6,*)'  [1] J. S. Jin, S. Welack, J. Y. Luo, X. Q. Li, P. Cui, R. X. Xu,       '
write(6,*)'      and Y. J. Yan, J. Chem. Phys. 126, 134113 (2007)                   '
write(6,*)'  [2] J. S. Jin, X. Zheng, and Y. J. Yan, J. Chem. Phys. 128,            '
write(6,*)'      234703 (2008)                                                      '
write(6,*)'  [3] X. Zheng, J. S. Jin, and Y. J. Yan, J. Chem. Phys. 129,            '
write(6,*)'      184112 (2008)                                                      '
write(6,*)'  [4] X. Zheng, J. S. Jin, and Y. J. Yan, New J. Phys. 10,               '
write(6,*)'      093016 (2008)                                                      '
write(6,*)'  [5] X. Zheng, J. Y. Luo, J. S. Jin, and Y. J. Yan, J. Chem. Phys.      '
write(6,*)'      130, 124508 (2009)                                                 '
write(6,*)'  [6] X. Zheng, J. S. Jin, S. Welack, M. Luo, and Y. J. Yan,             '
write(6,*)'      J. Chem. Phys. 130, 164708 (2009)                                  '
write(6,*)'  [7] J. Hu, R. X. Xu, and Y. J. Yan, J. Chem. Phys. 133, 101106 (2010)  '
write(6,*)'  [8] J. Hu, M. Luo, F. Jiang, R. X. Xu, and Y. J. Yan,                  '
write(6,*)'      J. Chem. Phys. 134, 244106 (2011)                                  '
write(6,*)'  [9] Z. Li, N. H. Tong, X. Zheng, J. H. Wei, J. Hu, and Y. J. Yan,      '
write(6,*)'      Phys. Rev. Lett. 109, 266403 (2012)                                '
write(6,*)' [10] S. Wang, X. Zheng, J. Jin, and Y. J. Yan,                          '
write(6,*)'      Phys. Rev. B 88, 035129 (2013)                                     '
write(6,*)' [11] D. Hou, R. Wang, X. Zheng, N. H. Tong, J. H. Wei, and Y. J. Yan,   '
write(6,*)'      Phys. Rev. B 90, 045141 (2014)                                     '
write(6,*)' [12] D. Hou, S. Wang, R. Wang, L. Ye, R. X. Xu, X. Zheng, and Y. J. Yan,'
write(6,*)'      J. Chem. Phys. 142, 104112 (2015)                                  '
write(6,*)' [13] L. Ye, X. Wang, D. Hou, R. X. Xu, X. Zheng, and Y. J. Yan,         '
write(6,*)'      WIREs Comput. Mol. Sci. 6, 608 (2016)                              '
write(6,*)' [14] H.-D. Zhang, Q. Qiao, R. X. Xu, X. Zheng, and Y. J. Yan,           '
write(6,*)'      J. Chem. Phys. 147, 044105 (2017)                                  '
write(6,*)' [15] L. Ye, H. D. Zhang, Y. Wang, X. Zheng, and Y. J. Yan,              '
write(6,*)'      J. Chem. Phys. 147, 074111 (2017)                                  '
write(6,*)' [16] L. Han, H.-D. Zhang, X. Zheng, and Y. J. Yan,                      '
write(6,*)'      J. Chem. Phys. 148, 234108 (2018)                                  '
write(6,*)' [17] L. Cui, H.-D. Zhang, X. Zheng, R.-X. Xu, and Y. J. Yan,            '
write(6,*)'      J. Chem. Phys. 151, 024110 (2019)                                  '
write(6,*)' [18] H.-D. Zhang, L. Cui, H. Gong, R.-X. Xu, X. Zheng, and Y. J. Yan,   '
write(6,*)'      J. Chem. Phys. 152, 064107 (2020)                                  '
write(6,*)' [19] D. Zhang, X. Ding, H.-D. Zhang, X. Zheng, and Y. J. Yan,           '
write(6,*)'      Chin. J. Chem. Phys. 34, 905 (2021)                                '
write(6,*)' [20] Z.-H. Chen, Y. Wang, X. Zheng, R.-X. Xu, and Y. J. Yan,            '
write(6,*)'      J. Chem. Phys. 156, 221102 (2022)                                  '
write(6,*)' [21] X. Ding, D. Zhang, L. Ye, X. Zheng, and Y. J. Yan,                 '
write(6,*)'      J. Chem. Phys. 157, 224107 (2022)                                  '
write(6,*)' [22] D. Zhang, L. Zuo, L. Ye, Z.-H. Chen, Y. Wang, R.-X. Xu, X. Zheng,  '
write(6,*)'      and Y. J. Yan, J. Chem. Phys. 158, 014106 (2023)                   '
write(6,*)
write(6,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
write(6,*)
write(6,1000)date(7:8),date(5:6),date(1:4),time(1:2),time(3:4),time(5:6)
!
1000 format(/,' Job started on ', /,                                     &
            ' ',A2,"/",A2,"/",A4,"  ",A2,":",A2,":",A2,/)
1001 format(15x, 'VERSION', 2x, f8.3)
!
write(6,*)
write(6,*)'genoutput: original input file ========================================> '
open(5, status='unknown')
rewind(5)
do 
  read(5, '(A120)', iostat=istat) line
  if (istat .eq. 0) then
     write(6, '(A120)') line
  else 
     exit
  end if
end do
write(6,*)'genoutput: <=========================================== end of input file'
call flush(6)
!
rewind(5)
!
end subroutine genoutput
