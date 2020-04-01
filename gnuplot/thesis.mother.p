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
set terminal postscript enhanced eps color font "Times-Roman,14"

mother="../diag/mother.signal.2867.dat"
field="../diag/fields.2867.dat"

rps="../diag/findmodes.right.18.dat"
lps="../diag/findmodes.left.19.dat"

lpw="../diag/power.left.dat"
rpw="../diag/power.right.dat"

# FIELD
set xlabel "x [{/Symbol l}_D]"
set ylabel "{/Symbol r}(x)"
set ylabel "E(x)"
set ylabel "E(x) [{/Symbol k} T/e{/Symbol l}_D]"

set key off
set size 1.0,0.5
set output "~/thesis/figures/mother.fullsignal.eps"
set title "Fullsize signal at t=92.15 {/Symbol w}_{pe}^{-1}"
plot field using 1:4 with lines 
unset size


# SIGNAL
set xlabel "x [{/Symbol l}_D]"
set ylabel "{/Symbol r}(x)"
set ylabel "E(x) [{/Symbol k} T/e{/Symbol l}_D]"

set xtics 50
set key top left box
set output "~/thesis/figures/mother.signal.eps"
set size 1,0.5
set multiplot layout 1,2
set size 0.5,0.5
set title "Left"
plot    mother using 1:3 with lines title "Original", \
        mother using 1:5 with lines title "Windowed"
set size 0.5,0.5
set title "Right"
plot    mother using 2:4 with lines title "Original", \
        mother using 2:6 with lines title "Windowed"
unset multiplot
unset size
set xtics auto

# POWER SPECTRUM
set xrange [0:0.3]
set ylabel "Power"
set xlabel "Wavenumber k [{/Symbol l}_D^{-1}]"

set key top right box
set output "~/thesis/figures/mother.ps.eps"
set multiplot layout 1,2 title "Accumulated power spectrum t=92.15 to t=102.35 {/Symbol w}_{pe}^{-1}"
set title "Left"
plot    lps using ($1*0.03351):2 with lines title "Power spectrum", \
        lps using ($1*0.03351):3 with lines title "Mean", \
        lps using ($1*0.03351):4 with lines title "Tolerance"

set title "Right"
plot    rps using ($1*0.03351):2 with lines title "Power spectrum", \
        rps using ($1*0.03351):3 with lines title "Mean", \
        rps using ($1*0.03351):4 with lines title "Tolerance"

unset multiplot
set autoscale

# TRACKING
set xlabel "Time [{/Symbol w}_{pe}^{-1}]"
set ylabel "Power"

set key left top box
set output "~/thesis/figures/mother.tracking.eps"

load 'results.power.p'

#set multiplot layout 1,2 title "Tracking of dominating modes"
#set title "Left"
#plot    lpw using ($1*0.05):2 with lines title "k=1", \
#        lpw using ($1*0.05):3 with lines title "k=2"
#
#set title "Right"
#plot    rpw using ($1*0.05):2 with lines title "k=1"
#unset multiplot







set out
set term wxt
unset title
unset label
unset xlabel
unset ylabel
unset size

