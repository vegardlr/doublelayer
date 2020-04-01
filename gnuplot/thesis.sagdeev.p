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
set terminal postscript enhanced eps color font "Times-Roman,24"
set output "~/thesis/figures/initDLsagdeev.eps"







L=250
data="../diag/supermethod.dat"






#set multiplot layout 3,1
set autoscale
# set ylabel "{/Symbol r}({/Symbol f}) & V({/Symbol f})" 
unset ylabel
set xlabel "{/Symbol f}"
set xtics ("0" 0, "{/Symbol f}_{DL}/2" 100,"{/Symbol f}_{DL}" 200)
set ytics (-1,0,1)
set yrange [-1.1:1.1]
set size 1,0.7
set key bottom right box

set title "Charge density, Sagdeev potential and position of potential"
plot    data using 2:($6/.000015) with lines title "{/Symbol r}({/Symbol f})/{/Symbol r}_{max}", \
        data using 2:($3/0.0005) with lines title "V({/Symbol f})/V_{max}", \
        data using 2:($4/L) with lines title "x({/Symbol f})/L"




set out
set output "~/thesis/figures/initDLphi.eps"

set size 1,0.5
set autoscale
set xtics ("0" 0, "L" L)
set ytics ("0" 0, "{/Symbol f}_{DL}/2" 100,"{/Symbol f}_{DL}" 200)
set yrange [-10:210]
set key off
set title "Potential profile {/Symbol f}(x)"
set xlabel "Position [{/Symbol l}_D]"
plot data using 1:5 with lines ls 1











set out
set term wxt
unset title
unset label
unset xlabel
unset ylabel
unset size

