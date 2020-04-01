

set terminal postscript enhanced color font "Times-Roman,12"
set output "../plot/findmodes.ps"


set title "Powerspectrum of {/Symbol r}(x)"
set ylabel "Power"
set xlabel "Wavenumber k "
m='../diag/findmodes.dat'
plot    m using 1:2 with lines title "Powerspectrum", \
        m using 1:3 with lines title "Mean value", \
        m using 1:4 with lines title "Tolerance"


set out
set term wxt
set autoscale
unset xlabel
unset ylabel
unset title

