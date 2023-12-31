#!/bin/csh
#
setenv WORK /home/xzheng/qdt-fermion/test
#
set icount = 1
set imax   = 21
#
if (-e $WORK/input.in)    \rm -f $WORK/input.in
if (-e $WORK/vin.tmp)     \rm -f $WORK/vin.tmp
if (-e $WORK/result.out)  \rm -f $WORK/result.out 
if (-e $WORK/output.log)  \rm -f $WORK/output.log
if (-e $WORK/variables.sav)   \rm -f $WORK/variables.sav
#
touch $WORK/result.out
touch $WORK/output.log
#
while ("$icount" <= "$imax")
  echo $icount > $WORK/vin.tmp
  $WORK/makeinput.x
  cat                 $WORK/input.in >  $WORK/out
  $WORK/qdt_trans.x < $WORK/input.in >> $WORK/out
  cat                 $WORK/out      >> $WORK/output.log
  \grep "IVOUTPUT"  < $WORK/out      >> $WORK/result.out
  @ icount++
end



