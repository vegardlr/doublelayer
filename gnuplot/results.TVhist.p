                        ##########
                        # TVhist #
                        ##########
ititle="Ions"
idata="`cat marks/TVhistp11`" 
etitle="Electrons"
edata="`cat marks/TVhiste11`"


              #########################
                #       TVhist          #
                #########################
set xlabel '(v_{in}^3/|v_{in}|)  [v_{th}^2]'
set ylabel '{/Symbol t}[{/Symbol w}_{pe}^{-1}]'
set logscale zcb
set multiplot layout 1,2 title tittel
set key off
set title etitle
set autoscale
plot  edata with image
set title ititle
set autoscale 
plot  idata with image
unset multiplot



