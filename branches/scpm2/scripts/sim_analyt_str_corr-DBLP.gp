set term postscript eps 30 enhanced color
set style fill solid 1.00 border
set encoding iso_8859_1
set output "sim_analyt_str_corr-DBLP.eps"
#set xlabel "{/Symbol s}(x10^3)"
set xlabel "{/Symbol s} x 10^3"
set ylabel "sim-{/Symbol e}_{exp} x 10^{-4}"
set y2label "max-{/Symbol e}_{exp} x 10^{-4}"
set xtic 2
set y2tic 50
set yrange [0:25]
set y2range [0:200]
set ytics nomirror
set key left top
set xrange [:11]
plot "<awk '{print $1/1000,$2*10000,$3*10000,$4*10000}' sim_analyt_str_corr-DBLP.dat" using 1:2:3 title "sim-{/Symbol e}_{exp}" with errorbars pt 2 lc 1 axes x1y1, "<awk '{print $1/1000,$2*10000,$3*10000,$4*10000}' sim_analyt_str_corr-DBLP.dat" using 1:4 title "max-{/Symbol e}_{exp}" with points pt 7 lc 3 axes x1y2
