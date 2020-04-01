


set multiplot layout 2,3 columnsfirst title "Double layer simulations (t=`cat ../diag/pres.time.dat` [{/Symbol w}_{pe}^{-1}])"

set key off
set title "Electrostatic potential"
set xlabel "x [{/Symbol l}_D]"
set ylabel "{/Symbol f}(x) [kT/e]"
set xrange [-256:384]
plot "../diag/pres.fields.dat" using 1:3 with lines
set autoscale

set title "Potential leap (temporal evolution)"
set xlabel "t [{/Symbol w}_{pe}^{-1}]"
set ylabel "{/Symbol f}_{DL}(t) [kT/e]"
set xrange [0:100]
plot "../diag/pres.phidl.dat" with lines 
set autoscale

set title "Phase Space Electrons"
set xlabel "x [{/Symbol l}_D]"
set ylabel "v [v_{th}^e]"
set xrange [-256:384]
set yrange [-30:30]
plot '../diag/pres.ps.xe-ve.dat' with image

set title "Phase Space Ions"
set xlabel "x [{/Symbol l}_D]"
set ylabel "v [v_{th}^e]"
set yrange [-1.75:1.75]
plot '../diag/pres.ps.xp-vp.dat' with image
set autoscale

set pm3d at bs
set view 60,60
set ylabel "x [{/Symbol l}_D]"
set xlabel "t [{/Symbol w}_{pe}^{-1}]"

set xrange [0:100]
set yrange [-256:384]
set title "Charge density"
splot "../diag/pres.carpet.dat" using 1:2:3 with pm3d
set title "Potential"
splot "../diag/pres.carpet.dat" using 1:2:4 with pm3d
set autoscale

unset pm3d
unset multiplot

