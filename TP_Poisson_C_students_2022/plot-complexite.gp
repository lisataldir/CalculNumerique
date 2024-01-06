set xlabel 'n'         
set ylabel 'performance (ns)' 

plot 'time0.dat' with title "dgbtrf + dgbtrs", \
    'time1.dat' with title " dgbtrftridiag + dgbtrs", \
    'time2.dat' with title "dgbsv"

pause -1