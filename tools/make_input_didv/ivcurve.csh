#!/bin/csh
#
setenv WORK /home/xzheng/WORK/HEOM/test_low-T
setenv BIN  /home/xzheng/WORK/HEOM/bin
#
set icount = 1 
set imax   = 33   
#
if (-e $WORK/input.in)    \rm -f $WORK/input.in
if (-e $WORK/vin.tmp)     \rm -f $WORK/vin.tmp
if (-e $WORK/result.out)  \rm -f $WORK/result.out 
if (-e $WORK/output.log)  \rm -f $WORK/output.log
#
touch $WORK/result.out
touch $WORK/output.log
#
while ("$icount" <= "$imax")
  echo $imax   >  $WORK/vin.tmp
  echo $icount >> $WORK/vin.tmp
  $WORK/makeinput.x
  $BIN/heom_quick.x < $WORK/input.in >  $WORK/out
  cat                 $WORK/out      >> $WORK/output.log
  \grep "IVOUTPUT"  < $WORK/out      >> $WORK/result.out
  @ icount++
end



