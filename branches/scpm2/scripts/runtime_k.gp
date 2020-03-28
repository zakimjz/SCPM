set term postscript eps 25 enhanced color
set style fill solid 1.00 border
set encoding iso_8859_1
set output "runtime_k.eps"
set xrange[0:17]
set yrange[:800000]
set xtic 2
set pointsize 2.5
set xlabel "k"
set ylabel "runtime (sec)"
set log y
plot 'runtime_k.dat' using 1:2 title "SCPM-DFS" with linespoints pt 7 lc 2 lt 1 lw 2,'runtime_k.dat' using 1:3 title "Naive" with linespoints pt 9 lc 3 lt 1 lw 2

