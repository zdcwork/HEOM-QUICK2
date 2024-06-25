#!/bin/csh
#
setenv WORK /home/jinjs/WORK/HEOM/test_ivtd
setenv BIN  /home/jinjs/WORK/HEOM/bin
#
set icount = 1 
set imax   = 10   
#
if (-e $WORK/input.in)    \rm -f $WORK/input.in
if (-e $WORK/vin.tmp)     \rm -f $WORK/vin.tmp
if (-e $WORK/result.out)  \rm -f $WORK/result.out 
if (-e $WORK/output.log)  \rm -f $WORK/output.log
if (-e $WORK/rho.out)     \rm -f $WORK/rho.out
#
touch $WORK/result.out
touch $WORK/output.log
touch $WORK/rho.out
#
while ("$icount" <= "$imax")
  echo $imax   >  $WORK/vin.tmp
  echo $icount >> $WORK/vin.tmp
  $WORK/makeinput.x
  $BIN/heom_fermion.x < $WORK/input.in >  $WORK/out
  cat                 $WORK/out      >> $WORK/output.log
  \grep " TDIV  "  < $WORK/out      >> $WORK/result.out
  \grep " TDRHO "  < $WORK/out      >> $WORK/rho.out   
  @ icount++
end



