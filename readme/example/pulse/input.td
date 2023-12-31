    1
    3
    2
    1
    2
    2
3.d0 3.d0
0.06d0  0.06d0 0.06d0 0.06d0
0.002d0 0.002d0
0.01d0   0.01d0    -0.01d0    -0.01d0
10.00
0.01

$para1  eup=-0.60  edown=-0.60  uu= 1.2d0    fixdot=F $end

$field  fieldtype= 2  $end
 1.0d-20  1.0d-20  1.0d-20   1.0d-20
pulsed_energy_shift  2  1
 1.d0     2.d0     1.0d0

$jobinfo  lsparse=T  psdfff=T  itype_fff=3  $end

$converge maxit0=30000 crit=1.d-7 $end

$adiabatic lad=.true. $end  

