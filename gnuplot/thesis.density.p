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

set terminal postscript enhanced eps color font "Times-Roman,10"
set output ".../plot/thesis.density.eps"






L=250
add=50

                #######################
                #       DENSITY       #
                #######################
p0="`cat marks/density0`"
plot0="p0 u 1:2 w l t 'Acc. electrons', p0 u 1:3 w l t 'Refl. electrons', p0 u 1:4 w l t'Accel. ions', p0 u 1:5 w l t 'Refl. ions'"
total0="p0 u 1:6 w l t 'Electrons', p0 u 1:7 w l t 'Ions'"


fields0="`cat marks/fields0`"
chargedensity = "fields0 using 1:2 with lines"
potential = "fields0 using 1:3 with lines"

unset xlabel
set xrange [-add:L+add]
set xtics ("0" 0, "L" L)
unset ytics

#set autoscale
#set ytics auto
#set xtics auto

set size 0.5,1
set multiplot layout 4,1 columnsfirst
set size 0.5,0.25
set key top center box title 'Population/Specie'
set ylabel "n_{sp}(x)"
set yrange [-100:5000]
plot @plot0
set ylabel "{/Symbol S}_p n_{sp}(x)"
set key top center box title 'Specie'
set size 0.5,0.25
set yrange [1000:5000]
plot @total0
set key off
set ylabel "{/Symbol r}(x)"
set size 0.5,0.25
set yrange [-0.035:0.035]
plot @chargedensity
set xlabel 'Position [{/Symbol l}_D]'
set ylabel "{/Symbol f}(x)"
set size 0.5,0.25
set yrange [0:220]
plot @potential

unset multiplot
unset size















set out
set term wxt
unset title
unset label
unset xlabel
unset ylabel


