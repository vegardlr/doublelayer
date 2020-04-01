                       ##############
                        # PHASESPACE #
                        ##############
psetitle1="Electrons t=`cat marks/time6` [{/Symbol w}_{pe}^{-1}]"
psedata1 ="`cat marks/psxe-ve6`"
psetitle2="Electrons t=`cat marks/time10` [{/Symbol w}_{pe}^{-1}]"
psedata2 ="`cat marks/psxe-ve10`"
psititle1="Ions t=`cat marks/time6` [{/Symbol w}_{pe}^{-1}]"
psidata1 ="`cat marks/psxp-vp6`"
psititle2="Ions t=`cat marks/time10` [{/Symbol w}_{pe}^{-1}]"
psidata2 ="`cat marks/psxp-vp10`"


set ylabel 'v [v_{th}]'
set xlabel 'x [{/Symbol l}_D]'
set multiplot layout 2,2 title tittel

set key off

set title psetitle1
plot psedata1 with image

set title psetitle2
plot psedata2 with image

set title psititle1
plot psidata1 with image

set title psititle2
plot psidata2 with image

unset multiplot

