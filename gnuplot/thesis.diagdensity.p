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

set terminal postscript enhanced eps color font "Times-Roman,12"
set output ".../plot/thesis.diagdensity.eps"






L=250
add=50

                #######################
                #       DENSITY       #
                #######################
p0="`cat marks/density0`"
plot0="p0 u 1:2 w l t 'Free e', p0 u 1:3 w l t 'Trapped e', p0 u 1:4 w l t 'Free i', p0 u 1:5 w l t 'Trapped i'"
total0="p0 u 1:6 w l t 'Electrons', p0 u 1:7 w l t 'Ions'"

set xlabel 'Position ({/Symbol l}_D)'
set xtics
set ytics

set key right center box
set size 1,0.5
set multiplot layout 1,2
set size 0.5,0.5
set title "Species"
plot @total0
set size 0.5,0.5
set title "Populations"
set yrange [-100:5000]
plot @plot0
unset multiplot


set key top center box title 'Population/Specie'
set ylabel "n_s^p(x)"
set yrange [0:500]
set ylabel "{/Symbol S}_p n_s^p(x)"
set key top center box title 'Specie'
set size 0.5,0.25
set yrange [100:500]
set key off
set ylabel "{/Symbol r}(x)"
set size 0.5,0.25
set yrange [-0.035:0.035]
set ylabel "{/Symbol f}(x)"
set size 0.5,0.25
set yrange [0:220]














set out
set term wxt
unset title
unset label
unset xlabel
unset ylabel


