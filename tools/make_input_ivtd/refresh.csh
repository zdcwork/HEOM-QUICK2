#!/bin/csh
ifort -o makeinput.x makeinput.f90
\cp -f *.x ivcurve.csh ../../work
