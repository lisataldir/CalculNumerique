set xlabel 'Itérations'
set ylabel 'Erreur'
set logscale y 10

plot 'resvec.dat' with linespoints title 'Erreur par itération'

pause -1