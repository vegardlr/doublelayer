
set terminal postscript portrait enhanced color font "Times-Roman,12"
set terminal postscript enhanced eps color font "Times-Roman,12"
set output "../plot/trace.ps"


unset grid
set key off
set title "v>0"
set xlabel "x"
set ylabel "v"
t="../diag/trace.dat"

set output "../plot/trace1.eps"
set multiplot layout 3,2 columnsfirst 
plot t u 2:3 w l, t u 4:5 w l, t u 6:7 w l, t u 8:9 w l, t u 10:11 w l, \
     t u 12:13 w l, t u 14:15 w l, t u 16:17 w l, t u 18:19 w l, \
     t u 20:21 w l 

set xtics ("0" 0, "102.375" 2049); set xrange [0:2049]
set xlabel "t"
set ylabel "x(t)"
plot t u 1:2 w l, t u 1:4 w l, t u 1:6 w l, t u 1:8 w l, t u 1:10 w l, \
        t u 1:12 w l,t u 1:14 w l,t u 1:16 w l,t u 1:18 w l,t u 1:20 w l 

set ylabel "v(t)"
plot t u 1:3 w l, t u 1:5 w l, t u 1:7 w l, t u 1:9 w l, t u 1:11 w l, \
        t u 1:13 w l,t u 1:15 w l,t u 1:17 w l,t u 1:19 w l,t u 1:21 w l 

set xtics; set autoscale
#unset multiplot
#set output "../plot/trace2.eps"
#set multiplot layout 3,1 columnsfirst 

set title "v<0"
set xlabel "x"
set ylabel "v"
plot t u 22:23 w l, t u 24:25 w l, t u 26:27 w l, t u 28:29 w l, \
        t u 30:31 w l, t u 32:33 w l, t u 34:35 w l, t u 36:37 w l, \
        t u 38:39 w l, t u 40:41 w l


set xtics ("0" 0, "102.375" 2049); set xrange [0:2049]
set xlabel "t"
set ylabel "x(t)"
plot t u 1:22 w l, t u 1:24 w l, t u 1:26 w l, t u 1:28 w l, t u 1:30 w l, \
        t u 1:32 w l,t u 1:34 w l,t u 1:36 w l,t u 1:38 w l,t u 1:40 w l 

set ylabel "v(t)"
plot t u 1:23 w l, t u 1:25 w l, t u 1:27 w l, t u 1:29 w l, t u 1:31 w l, \
        t u 1:33 w l,t u 1:35 w l,t u 1:37 w l,t u 1:39 w l,t u 1:41 w l 

set xtics; set autoscale
unset multiplot

set output "../plot/trace3.eps"
#set multiplot layout 2,1
#
#set pm3d at bs
#set view 70,160
#set grid
#set yrange [-100:-80]
#set ylabel "x"
#set xlabel "t"
#set title "phi(x)"
#splot "../diag/carpet.dat" using ($1/0.05):2:4 with pm3d

set xtics ("0" 150, "102.375" 2049); set xrange [0:2049]
set autoscale
set yrange [0:2]
set xrange [105:2049]
set pm3d at bs
set title "Influx histogram of reflected ions"
set xlabel "Time"
set ylabel "Velocity"
#set xrange [100:500]
#set yrange [0:3]
set view 70,70
splot "../diag/influxhist.dat" with pm3d

unset pm3d
unset multiplot

set out
set term wxt
set autoscale
unset xlabel
unset ylabel
unset title

