set macros
set grid
set terminal postscript portrait enhanced color
set output "../plot/field.ps"

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

options1="using 1:2 with lines title"
options2="using 1:3 with lines title"
set multiplot layout 2,1
plot    fields0 @options1 @time0
plot    fields0 @options2 @time0
unset multiplot
set multiplot layout 2,1
plot    fields1 @options1 @time1
plot    fields1 @options2 @time1
unset multiplot
set multiplot layout 2,1
plot    fields2 @options1 @time2
plot    fields2 @options2 @time2
unset multiplot
set multiplot layout 2,1
plot    fields3 @options1 @time3
plot    fields3 @options2 @time3
unset multiplot
set multiplot layout 2,1
plot    fields4 @options1 @time4
plot    fields4 @options2 @time4
unset multiplot
set multiplot layout 2,1
plot    fields5 @options1 @time5
plot    fields5 @options2 @time5
unset multiplot
set multiplot layout 2,1
plot    fields6 @options1 @time6
plot    fields6 @options2 @time6
unset multiplot
set multiplot layout 2,1
plot    fields7 @options1 @time7
plot    fields7 @options2 @time7
unset multiplot
set multiplot layout 2,1
plot    fields8 @options1 @time8
plot    fields8 @options2 @time8
unset multiplot
set multiplot layout 2,1
plot    fields9 @options1 @time9
plot    fields9 @options2 @time9
unset multiplot
set multiplot layout 2,1
plot    fields10 @options1 @time10
plot    fields10 @options2 @time10
unset multiplot

set out
set terminal wxt
