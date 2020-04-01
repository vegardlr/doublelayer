                        ###########
                        # DENSITY #
                        ###########
p0="`cat marks/density0`"
p1="`cat marks/density2`"
title1="t=`cat marks/time0` vs t=`cat marks/time2` [{/Symbol w}_{pe}^{-1}]"
p2="`cat marks/density4`"
title2="t=`cat marks/time0` vs t=`cat marks/time4` [{/Symbol w}_{pe}^{-1}]"
p3="`cat marks/density6`"
title3="t=`cat marks/time0` vs t=`cat marks/time6` [{/Symbol w}_{pe}^{-1}]"
 

plot0="p0 u 1:2 w l t 'ea', p0 u 1:3 w l t 'er', p0 u 1:4 w l t'ia', p0 u 1:5 w l t 'ir'"
plot1="p1 u 1:2 w l t 'ea', p1 u 1:3 w l t 'er', p1 u 1:4 w l t'ia', p1 u 1:5 w l t 'ir'"
plot2="p2 u 1:2 w l t 'ea', p2 u 1:3 w l t 'er', p2 u 1:4 w l t'ia', p2 u 1:5 w l t 'ir'"
plot3="p3 u 1:2 w l t 'ea', p3 u 1:3 w l t 'er', p3 u 1:4 w l t'ia', p3 u 1:5 w l t 'ir'"
total0="p0 u 1:6 w l t 'e ', p0 u 1:7 w l t 'i '"
total1="p1 u 1:6 w l t 'e ', p1 u 1:7 w l t 'i '"
total2="p2 u 1:6 w l t 'e ', p2 u 1:7 w l t 'i '"
total3="p3 u 1:6 w l t 'e ', p3 u 1:7 w l t 'i '"

set autoscale
set xlabel 'x [{/Symbol l}_D]'
set multiplot layout 2,3 columnsfirst title tittel

set key outside left center box title 'Specie' samplen 2
set title title1
#set yrange[-100:5500]
set ylabel 'n_{sp}(x) [n_0]'
plot @plot0, @plot1
set ylabel 'n_{s}(x) [n_0]'
#set yrange [1000:5500]
plot @total0, @total1

set key off
set title title2
#set yrange[-100:5500]
set ylabel 'n_{sp}(x) [n_0]'
plot @plot0, @plot2
#set yrange [1000:5500]
set ylabel 'n_{s}(x) [n_0]'
plot @total0, @total2

set title title3
#set yrange[-100:5500]
set ylabel 'n_{sp}(x) [n_0]'
plot @plot0, @plot3
#set yrange [1000:5500]
set ylabel 'n_{s}(x) [n_0]'
plot @total0, @total3


unset multiplot
set autoscale

