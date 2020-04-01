HODE = "`pwd | tail -c 13`"

TITTEL= sprintf(HODE,FIL)
set nokey
set yrange [-10:25]
set xrange [-128:256]

 FIL="vND/diag/ps.xe-ve.0.dat" 
 set title FIL
 plot FIL with image
 !sleep 1

 FIL="vND/diag/ps.xe-ve.50.dat" 
 set title FIL
 plot FIL with image
 !sleep 1

 FIL="vND/diag/ps.xe-ve.100.dat" 
 set title FIL
 plot FIL with image
 !sleep 1

 FIL="vND/diag/ps.xe-ve.150.dat" 
 set title FIL
 plot FIL with image
 !sleep 1

 FIL="vND/diag/ps.xe-ve.200.dat" 
 set title FIL
 plot FIL with image
 !sleep 1

 FIL="vND/diag/ps.xe-ve.250.dat" 
 set title FIL
 plot FIL with image
 !sleep 1

 FIL="vND/diag/ps.xe-ve.300.dat" 
 set title FIL
 plot FIL with image
 !sleep 1

 FIL="vND/diag/ps.xe-ve.350.dat" 
 set title FIL
 plot FIL with image
 !sleep 1

 FIL="vND/diag/ps.xe-ve.400.dat" 
 set title FIL
 plot FIL with image
 !sleep 1

 FIL="vND/diag/ps.xe-ve.450.dat" 
 set title FIL
 plot FIL with image
 !sleep 1

 FIL="vND/diag/ps.xe-ve.500.dat" 
 set title FIL
 plot FIL with image
 !sleep 1

