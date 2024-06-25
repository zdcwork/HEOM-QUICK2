    1
    3
    3
    2
    2
    2
 0.500000E+001    0.500000E+001
 0.100000E-001    0.100000E-001    0.600000E-001    0.600000E-001
 0.640000E-003    0.640000E-003
 0.000000E+000    0.000000E+000    0.200000E-002    0.200000E-002
 0.200000E+004
 0.500000E-002

$para4 e1up=  -2.5000  e1down=  -2.5000  e2up=  -2.5000  e2down=  -2.5000
 uu1=   5.0000  uu2=   5.0000  u12=   0.0000  t12=   0.0000 J12=   0.0000 $end

$converge  maxit0= 20000   crit=  1.00E-007  $end

$jobinfo lsparse=.true.  psfjob=T  itype_psf=1  $end

$field  fieldtype= 0  $end
 0.100000E-019    0.100000E-019    0.100000E-019    0.100000E-019

$resume  lresume=      T  nresume=   20000  icont=  0  $end

$spin3d  lspin3d=.true.  $end

$zfs  lzfs=.true.  d_xx= -0.0040  d_yy = -0.0070  d_zz = -0.0110  $end

