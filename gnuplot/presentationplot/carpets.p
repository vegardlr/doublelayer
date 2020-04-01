set view 60,60
set pm3d at bs
set autoscale
set xlabel 'Time ({/Symbol w}_{pe}^{-1})'
set ylabel 'Position ({/Symbol l}_D)'
set key off

set terminal postscript enhanced color
set output '../plot/carpets.ps'


set title "Charge density\n `cat ../diag/params.dat`"
set zlabel "{/Symbol r}(x)"
splot '../diag/carpet.dat' using 1:2:3 with pm3d

set title "Potential\n `cat ../diag/params.dat`"
set zlabel "{/Symbol f}(x)"
splot '../diag/carpet.dat' using 1:2:4 with pm3d

set title "Electric field\n `cat ../diag/params.dat`"
set zlabel "E(x)"
splot '../diag/carpet.dat' using 1:2:5 with pm3d

set out
set terminal wxt
unset xlabel
unset ylabel
unset zlabel
unset title
unset pm3d
