                        #########
                        # VDF   #
                        #########

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


vhe0="`cat marks/vdfe0`"
vhe1="`cat marks/vdfe1`"
vhe2="`cat marks/vdfe2`"
vhe3="`cat marks/vdfe3`"
vhe4="`cat marks/vdfe4`"
vhe5="`cat marks/vdfe5`"
vhe6="`cat marks/vdfe6`"
vhe7="`cat marks/vdfe7`"
vhe8="`cat marks/vdfe8`"
vhe9="`cat marks/vdfe9`"
vhe10="`cat marks/vdfe10`"

vhp0="`cat marks/vdfp0`"
vhp1="`cat marks/vdfp1`"
vhp2="`cat marks/vdfp2`"
vhp3="`cat marks/vdfp3`"
vhp4="`cat marks/vdfp4`"
vhp5="`cat marks/vdfp5`"
vhp6="`cat marks/vdfp6`"
vhp7="`cat marks/vdfp7`"
vhp8="`cat marks/vdfp8`"
vhp9="`cat marks/vdfp9`"
vhp10="`cat marks/vdfp10`"



                #########################       
                #       VDF             #
                #########################
set xlabel 'v [v_{th}]'
set multiplot layout 2,2 rowsfirst title tittel


set key off
options = "using 1:2 with histeps title"
set title "Electrons x=lx"
plot    vhe0 @options @time0, vhe1 @options @time1, vhe2 @options @time2, vhe3 @options @time3,         vhe4 @options @time4, vhe5 @options @time5,         vhe6 @options @time6, vhe7 @options @time7,         vhe8 @options @time8, vhe9 @options @time9,         vhe10 @options @time10



set key outside left center box samplen 2 width 0 title 't [{/Symbol w}_{pe}^{-1}]' 
options = "using 1:5 with histeps title"
set title "Electrons x=ux"
plot    vhe0 @options @time0, vhe1 @options @time1,         vhe2 @options @time2, vhe3 @options @time3,         vhe4 @options @time4, vhe5 @options @time5,         vhe6 @options @time6, vhe7 @options @time7,         vhe8 @options @time8, vhe9 @options @time9,         vhe10 @options @time10



set key off
options = "using 1:2 with histeps title"
set title "Ions x=lx"
plot    vhp0 @options @time0, vhp1 @options @time1,         vhp2 @options @time2, vhp3 @options @time3,         vhp4 @options @time4, vhp5 @options @time5,         vhp6 @options @time6, vhp7 @options @time7,         vhp8 @options @time8, vhp9 @options @time9,         vhp10 @options @time10



set key off
options = "using 1:5 with histeps title"
set title "Ions x=ux"
plot    vhp0 @options @time0, vhp1 @options @time1,         vhp2 @options @time2, vhp3 @options @time3,         vhp4 @options @time4, vhp5 @options @time5,         vhp6 @options @time6, vhp7 @options @time7,         vhp8 @options @time8, vhp9 @options @time9,         vhp10 @options @time10

unset multiplot




