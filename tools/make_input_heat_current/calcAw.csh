#!/bin/csh
#
setenv WORK ./
setenv BIN  /home/xzheng/WORK/HEOM/bin
#
set icount = 1 
set imax   = 26
#
if (-e $WORK/input.in)    \rm -f $WORK/input.in
if (-e $WORK/vin.tmp)     \rm -f $WORK/vin.tmp
if (-e $WORK/Iw.dat)  \rm -f $WORK/Iw.dat 
if (-e $WORK/JHw.dat)  \rm -f $WORK/JHw.dat 
if (-e $WORK/JEw.dat)  \rm -f $WORK/JEw.dat 
if (-e $WORK/output.log)  \rm -f $WORK/output.log
#
touch $WORK/Iw.dat
touch $WORK/JHw.dat
touch $WORK/JEw.dat
touch $WORK/output.log
#
while ("$icount" <= "$imax")
  echo $imax   >  $WORK/vin.tmp
  echo $icount >> $WORK/vin.tmp
  $WORK/makeinput.x
  $BIN/heom_fermion.x < $WORK/input.in >  $WORK/out
  cat                 $WORK/out      >> $WORK/output.log
  \grep "IW"        < $WORK/out      >> $WORK/Iw.dat
  \grep "JHW"       < $WORK/out      >> $WORK/JHw.dat
  \grep "JEW"       < $WORK/out      >> $WORK/JEw.dat
  @ icount++
end



