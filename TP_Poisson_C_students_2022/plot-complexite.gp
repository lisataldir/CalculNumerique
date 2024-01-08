set xlabel 'n'
set ylabel 'performance (ms)'

plot 'time0.dat' with lines smooth cspline title "dgbtrf + dgbtrs", \
    'time1.dat' with lines smooth cspline title "dgbtrftridiag + dgbtrs sans tmp", \
    'time1bis.dat' with lines smooth cspline title "dgbtrftridiag + dgbtrs avec tmp", \
    'time2.dat' with lines smooth cspline title "dgbsv"

pause -1