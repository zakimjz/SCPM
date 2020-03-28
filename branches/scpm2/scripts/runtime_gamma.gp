set term postscript eps 25 enhanced color
set style fill solid 1.00 border
set encoding iso_8859_1
set output "runtime_gamma.eps"
set xrange[0.45:1.05]
set yrange[:800000]
set xtic 0.1
set pointsize 2.5
set xlabel "{/Symbol g}_{min}"
set ylabel "runtime (sec)"
set log y
plot 'runtime_gamma.dat' using 1:2 title "SCPM-BFS" with linespoints pt 5 lc 1 lt 1 lw 2,'runtime_gamma.dat' using 1:4 title "Naive" with linespoints pt 9 lc 3 lt 1 lw 2,'runtime_gamma.dat' using 1:3 title "SCPM-DFS" with linespoints pt 7 lc 2 lt 1 lw 2

