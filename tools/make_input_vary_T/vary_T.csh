#!/bin/csh
#
setenv WORK ./
setenv BIN  /home/xzheng/WORK/HEOM/bin
#
set icount = 1 
set imax   = 21
#
if (-e $WORK/input.in)    \rm -f $WORK/input.in
if (-e $WORK/vin.tmp)     \rm -f $WORK/vin.tmp
if (-e $WORK/output.log)  \rm -f $WORK/output.log
if (-e $WORK/jup.dat)     \rm -f $WORK/jup.dat
if (-e $WORK/jdn.dat)     \rm -f $WORK/jdn.dat
if (-e $WORK/sz.dat)      \rm -f $WORK/sz.dat
if (-e $WORK/jdiff.dat)   \rm -f $WORK/jdiff.dat
#
touch $WORK/output.log
touch $WORK/jup.dat
touch $WORK/jdn.dat
touch $WORK/sz.dat
#
while ("$icount" <= "$imax")
  echo $imax   >  $WORK/vin.tmp
  echo $icount >> $WORK/vin.tmp
  $WORK/makeinput.x
  $BIN/heom_fermion.x < $WORK/input.in >  $WORK/out
  cat                 $WORK/out      >> $WORK/output.log
  \grep 'j_up  (ialf =            2  )=' $WORK/out >> $WORK/jup.dat
  \grep 'j_down(ialf =            2  )=' $WORK/out >> $WORK/jdn.dat
  \grep '< Sz >' $WORK/out >> $WORK/sz.dat
  @ icount++
end
#
paste $WORK/jup.dat $WORK/jdn.dat | awk '{print $12 - $6}' > $WORK/jdiff.dat
# 
exit 0



