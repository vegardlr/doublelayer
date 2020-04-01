                        ##################
                        # POWER TRACKING #
                        ##################

unset multiplot


set ylabel "Power"
set xlabel "t [{/Symbol w}_{pe}^{-1}]"
set logscale y
set key inside left top box title "k"




#set size 1,0.5
#set multiplot layout 1,2

set title "Left"
#set size 0.5,0.5
load '../diag/power.left.p'

#set size 0.5,0.5
#set title "Right"
#load '../diag/power.right.p'

unset multiplot

unset size
unset logscale y


