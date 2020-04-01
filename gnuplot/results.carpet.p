set view 60,60
set pm3d at bs
set autoscale
set xlabel 't [{/Symbol w}_{pe}^{-1}]'
set ylabel 'x [{/Symbol l}_D]'
set cbtics
set key off


set title tittel
set zlabel "{/Symbol r}(x,t)           \n [en_0]           "
splot '../diag/carpet.dat' using 1:2:3 with pm3d

#set title "Potential\n `cat ../diag/params.dat`"
#set zlabel "{/Symbol f}(x,t) [{/Symbol k}T/e]"
#splot '../diag/carpet.dat' using 1:2:4 with pm3d

#set title "Electric field\n `cat ../diag/params.dat`"
#set zlabel "E(x,t) [{/Symbol k}T/e{/Symbol l}_D]"
#splot '../diag/carpet.dat' using 1:2:5 with pm3d

unset xlabel
unset ylabel
unset zlabel
unset title
unset pm3d
