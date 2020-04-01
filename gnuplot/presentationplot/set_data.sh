#!/bin/sh
# Copies data to a filname plotted by Gnuplotscript runpres.p and
# presentatation.p. $1 is time of data

# Time is converted to physical units, given that dt=0.1
echo $((${1}/10)) > ../diag/pres.time.dat

# The copying
cp ../diag/fields_$1.dat ../diag/pres.fields.dat
cp ../diag/phidl.$1.dat ../diag/pres.phidl.dat
cp ../diag/ps.xe-ve.$1.dat ../diag/pres.ps.xe-ve.dat
cp ../diag/ps.xp-vp.$1.dat ../diag/pres.ps.xp-vp.dat
cp ../diag/carpet.$1.dat ../diag/pres.carpet.dat
