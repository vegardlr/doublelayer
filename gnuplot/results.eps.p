load 'script.header.p'





plotinitial     =1
plotfields      =1
plotcarpet      =1
plotdensity     =1
plotflux        =1
plotvdf         =1
plotps          =1
plotTVhist      =1
plotpower       =1




set terminal postscript enhanced eps color font "Times-Roman,11"

tittel="`cat ../diag/partnum.dat``cat ../diag/params.dat`"
carpettittel="`cat ../diag/partnum.dat`\n`cat ../diag/params.dat`"

if(plotinitial) set output "../plot/results.initial.eps"
if(plotinitial) load 'results.initial.p'
if(plotinitial) set out

if(plotfields)  set output "../plot/results.fields.eps"
if(plotfields)  load 'results.fields.p' 
if(plotfields)  set out

if(plotcarpet)  set output "../plot/results.carpet.eps"
if(plotcarpet)  load 'results.carpet.p'
if(plotcarpet)  set out

if(plotdensity) set output "../plot/results.density.eps"
if(plotdensity) load 'results.density.p'
if(plotdensity) set out

if(plotflux)    set output "../plot/results.flux.eps"
if(plotflux)    load 'results.flux.p'
if(plotflux)    set out

if(plotvdf)     set output "../plot/results.vdf.eps"
if(plotvdf)     load 'results.vdf.p'
if(plotvdf)     set out

if(plotps)      set output "../plot/results.phasespace.eps"
if(plotps)      load 'results.phasespace.p'
if(plotps)      set out

if(plotTVhist)  set output "../plot/results.TVhist.eps"
if(plotTVhist)  load 'results.TVhist.p'
if(plotTVhist)  set out

if(plotpower)   set output "../plot/results.power.eps"
if(plotpower)   load 'results.power.p'
if(plotpower)   set out


set out

load 'script.footer.p'

