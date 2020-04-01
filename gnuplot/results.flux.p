
                #################
                #       FLUX    #
                #################
doplot=0
if(doplot)        set multiplot layout 4,1 title "Flux \n`cat ../diag/partnum.dat`  `cat ../diag/params.dat`\n `cat ../diag/date.dat`"

if(doplot)        set xrange [`cat marks/time0`:`cat marks/time10`]
if(doplot)        flux="../diag/flux.dat"
if(doplot)        set xlabel "Time [{/Symbol w}_{pe}^{-1}]"
if(doplot)        set key bottom left box
if(doplot)        set title "Number of particles"
if(doplot)        plot    flux u 1:2 w l t "Electrons", flux u 1:3 w l t "Ions"
if(doplot)        set key top left box
if(doplot)        set title "Total flux, electrons"
if(doplot)        plot    flux u 1:4 w l t "In", flux u 1:5 w l t "Out"
if(doplot)        set title "Total flux, ions"
if(doplot)        plot    flux u 1:6 w l t "In", flux u 1:7 w l t "Out"
if(doplot)        set key off
if(doplot)        set title "Net charge change"
if(doplot)        plot    flux u 1:8 w l
if(doplot)        unset multiplot

f='../diag/flux2.dat'
set xlabel "t [{/Symbol w}_{pe}^{-1}]"
unset ylabel

set multiplot layout 2,1 title tittel
set key center left box width -3
set title 'Electrons'
plot    f u 1:2 w l t 'In acc',          \
        f u 1:3 w l t 'Out acc x=0',       \
        f u 1:4 w l t 'Out acc x=L',       \
        f u 1:5 w l t 'In refl',           \
        f u 1:6 w l t 'Out refl x=0',        \
        f u 1:7 w l t 'Out refl x=L'
set key center left box width -3
set title 'Ions'
plot    f u 1:8 w l t 'In acc',          \
        f u 1:9 w l t 'Out acc x=0',       \
        f u 1:10 w l t 'Out acc x=L',       \
        f u 1:11 w l t 'In refl',          \
        f u 1:12 w l t 'Out refl x=0',       \
        f u 1:13 w l t 'Out refl x=L'
unset multiplot
set autoscale


