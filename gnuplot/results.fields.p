                #########################       
                #       FIELDS          #
                #########################


time0="'`cat marks/time0`'"
time1="'`cat marks/time1`'"
time2="'`cat marks/time2`'"
time3="'`cat marks/time3`'"
time4="'`cat marks/time4`'"
time5="'`cat marks/time5`'"
time6="'`cat marks/time6`'"
time7="'`cat marks/time7`'"
time8="'`cat marks/time8`'"
time9="'`cat marks/time9`'"
time10="'`cat marks/time10`'"

fields0="`cat marks/fields0`"
fields1="`cat marks/fields1`"
fields2="`cat marks/fields2`"
fields3="`cat marks/fields3`"
fields4="`cat marks/fields4`"
fields5="`cat marks/fields5`"
fields6="`cat marks/fields6`"
fields7="`cat marks/fields7`"
fields8="`cat marks/fields8`"
fields9="`cat marks/fields9`"
fields10="`cat marks/fields10`"








set xlabel 'x [{/Symbol l}_D]'

set multiplot layout 2,2 title tittel

set key off
options="using 1:2 w l t"
set title "Charge density"
set ylabel '{/Symbol r}(x) [en_0]'
plot fields0 @options @time0, fields1 @options @time1, fields2 @options @time2, fields3 @options @time3, fields4 @options @time4, fields5 @options @time5, fields6 @options @time6, fields7 @options @time7, fields8 @options @time8, fields9 @options @time9, fields10 @options @time10


set key box samplen 2 width 0 title 't [{/Symbol w}_{pe}^{-1}]' outside left center
options="using 1:3 w l t"
set title "Potential"
set ylabel '{/Symbol f}(x) [{/Symbol k}T/e]'
plot    fields0 @options @time0, fields1 @options @time1, fields2 @options @time2, fields3 @options @time3, fields4 @options @time4, fields5 @options @time5, fields6 @options @time6, fields7 @options @time7, fields8 @options @time8, fields9 @options @time9, fields10 @options @time10

set key off
options="using 1:4 w l t"
set title "Electric field"
set ylabel 'E(x) [{/Symbol k}T/e{/Symbol l}_D]'
plot    fields0 @options @time0, fields1 @options @time1, fields2 @options @time2, fields3 @options @time3, fields4 @options @time4, fields5 @options @time5, fields6 @options @time6, fields7 @options @time7, fields8 @options @time8, fields9 @options @time9, fields10 @options @time10



set key off
set title "Evolution of {/Symbol f}_{DL} in time"
set xlabel "t [{/Symbol w}_{pe}^{-1}]"
set ylabel '{/Symbol f}_{DL}(t) [{/Symbol k}T/e]'
plot    '../diag/phidl.dat' w l t "phiDL(time)"

set key default
unset multiplot

