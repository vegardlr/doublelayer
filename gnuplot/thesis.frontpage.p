load 'script.header.p'

set terminal wxt enhanced
set terminal postscript enhanced eps color font "Times-Roman,11"
set output ".../plot/thesis.frontpage.eps"

set key off


set ylabel 'v [v_{the}]'
set xlabel 'x [{/Symbol l}_D]'
set xtics 150

psetitle2="Phase space diagram, electrons"
psedata2 ="`cat marks/psxe-ve10`"
psititle2="Phase space diagram, ions"
psidata2 ="`cat marks/psxp-vp10`"

set size 0.5,1
set multiplot layout 2,1

set size 0.5,0.5
set title psetitle2
plot psedata2 with image

set size 0.5,0.5
set title psititle2
plot psidata2 with image

unset multiplot
unset size

load 'script.footer.p'
