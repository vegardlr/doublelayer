
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
set output "~/thesis/figures/initDLspace.eps"

L=250


data="../diag/initialphi.dat"
set size 1,0.5
set autoscale
set xtics ("x_{lb}"  -L, "0" 0, "L" L, "x_{ub}" 2*L)
set ytics ("0" 0, "{/Symbol f}_{DL}/2" 100,"{/Symbol f}_{DL}" 200)
set yrange [-10:210]
set key off
set ylabel "{/Symbol f}(x)  [{/Symbol k}T/e]"
set xlabel "Position [{/Symbol l}_D]"
plot data with lines


set out
