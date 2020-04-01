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
set terminal postscript enhanced eps color font "Times-Roman,18"
set output "~/thesis/figures/printf.eps"
## set output "~/thesis/figures/printf_phi.eps"





   ###########################################################
   #               PLOTTING SECTION                          #
   ###########################################################
                #########################
                #       INITIAL         #
                #########################
L=250
PHI=200
set key off
set ytics auto
set xtics auto
unset xtics
unset ytics
unset ztics
#unset cbtics
set xtics out ("x_{lb}" -L,  "0" 0, "L" L, "x_{ub}" (L+L-2))
## set xtics ("0" 0, "{/Symbol f}_{DL}" PHI)
unset size
unset y2tics
unset y2label

set cbtics ("qqq" 0,"rrr" 1)

set size 1.0,0.6
set multiplot layout 1,2
#set pm3d at bs
set view 70,70
set ylabel 'Velocity'
set xlabel 'Position [{/Symbol l}_D]'
## set xlabel "Potential [{/symbol k}T/e]"


set ytics out ("-10 v_{the}" -10,  "0" 0, "10 v_{the}" 10)
set title "f_e(x,v)"
## set title "f_e({/Symbol f},v)"
set size 0.50,0.6
plot '../diag/printf.e.dat' with image

unset ylabel
set ytics ("-10 v_{thi}" -1,  "0" 0, "10 v_{thi}" 1)
set title "f_i(x,v)"
## set title "f_i({/Symbol f},v)"
set size 0.50,0.6
plot '../diag/printf.i.dat' with image
unset ylabel
unset pm3d
unset multiplot









set out
set term wxt
unset title
unset label
unset xlabel
unset ylabel


