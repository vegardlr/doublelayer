               #########################
                #       INITIAL         #
                #########################

fields0="`cat marks/fields0`"

s='../diag/supermethod.dat'
set xlabel 'x [{/Symbol l}_D]'


set multiplot layout 2,2 title tittel
#"Initial state\n`cat ../diag/partnum.dat` `cat ../diag/params.dat`\n `cat ../diag/date.dat`"

set key top left box
set title "{/Symbol r}(x)"
set ylabel "{/Symbol r}(x) [en_0]"
plot    s u 1:7 w l t 'expected', \
       fields0 u 1:2 w l t 'initialized'

set title "{/Symbol f}(x)"
set ylabel "{/Symbol f}(x) [{/Symbol k}T/e]"
plot    s u 1:5 w l t 'expected', \
        fields0 u 1:3 w l t 'initialized'
set key off

set ylabel 'v [v_{th}]'

set title "f_e(x,v)"
plot '../diag/printf.e.dat' with image

set title "f_i(x,v)"
plot '../diag/printf.i.dat' with image



unset ylabel
unset multiplot


