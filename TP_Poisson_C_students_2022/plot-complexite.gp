set xlabel 'n'
set ylabel 'performance (ms)'

plot 'time0.dat' with lines title "dgbtrf + dgbtrs", \
    'time1.dat' with lines title "dgbtrftridiag + dgbtrs", \
    'time2.dat' with lines title "dgbsv"

pause -1