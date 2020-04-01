# GNUPLOT file for plotting results based on fields from Double Layer
# simulations as a part of my masters thesis. -Vegard L. Rekaa
# Based on GNUPLOT 4.2

set macros                             # enable macros
set autoscale                          # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
unset title
unset xlabel
unset ylabel
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
unset multiplot

set terminal wxt enhanced
#set terminal postscript enhanced eps color font "Times-Roman,10"
set output ".../plot/thesis.densitytanh.eps"







L=100

                #######################
                #       DENSITY       #
                #######################


unset xlabel
set xtics ("0" 0, "L" L)
unset ytics
#set multiplot layout 4,1 columnsfirst title "Densities of a DL"
set key top center box title 'Population/Specie'
data="../diag/tanhphi.dat"
set ylabel "n_s^p(x)"
plot    data using 1:3 with lines title "Accel. e", \
        data using 1:4 with lines title "Refl. e", \
        data using 1:5 with lines title "Decel. e", \
        data using 1:6 with lines title "Accel. i", \
        data using 1:7 with lines title "Refl. i", \
        data using 1:8 with lines title "Decel. i"
!sleep 5
set ylabel "{/Symbol S}_p n_s^p(x)"
plot    data using 1:($3+$4+$5) with lines title "Electrons", \
        data using 1:($6+$7+$8) with lines title "Ions"
!sleep 5
set key off
set ylabel "{/Symbol f}(x)"
plot    data using 1:2 with lines title "Phi"
set xlabel 'Position ({/Symbol l}_D)'

#unset multiplot














set out
set term wxt
unset title
unset label
unset xlabel
unset ylabel


