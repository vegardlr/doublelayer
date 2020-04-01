load 'script.header.p'

plotinitial     =1
plotfields      =1
plotcarpet      =0
plotdensity     =1
plotflux        =1
plotvdf         =1
plotps          =1
plotTVhist      =1
plotpower       =1

set terminal postscript enhanced color font "Times-Roman,8"
set output "../plot/results.ps"

tittel="`cat ../diag/partnum.dat` `cat ../diag/params.dat`\n `cat ../diag/date.dat`"
carpettittel="`cat ../diag/partnum.dat`\n`cat ../diag/params.dat`"


if(plotinitial) load 'results.initial.p'

if(plotfields)  load 'results.fields.p' 

if(plotcarpet)  load 'results.carpet.p'

if(plotdensity) load 'results.density.p'

if(plotflux)    load 'results.flux.p'

if(plotvdf)     load 'results.vdf.p'

if(plotps)      load 'results.phasespace.p'

if(plotTVhist)  load 'results.TVhist.p'

if(plotpower)   load 'results.power.p'

set out

load 'script.footer.p'

