    1
    4
    2
    1
    2
    1
5.d0
0.24d0  0.24d0
0.0002d0
0.d0   0.d0
400.00d0
0.005d0

$para1  eup=-0.48d0  edown=-0.48d0  uu= 0.96d0    fixdot=F $end

$field  fieldtype= 0  $end
0.10000000E-19   0.10000000E-19

$jobinfo  lsparse=T  psfjob=T  itype_psf=33  $end

$converge maxit0=20000 crit=1.d-7 $end

$resume  icont=0  lresume=T   nresume=10000  $end

$adiabatic lad=.true.  $end

