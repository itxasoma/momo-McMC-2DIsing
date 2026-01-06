set terminal pdfcairo enhanced color

set output "jackerrs_Xi-C_vs_block.pdf"
data = "jackerr_Xi-C_vs_block.out"

set xrange [2:]
set logscale xy
set grid
set key bottom right 
set xlabel "block size"
set ylabel "{/Symbol s}"

plot data index 1 u 1:2 w linespoints pt 7 t "M", data index 0 u 1:3 w linespoints pt 7 t "E",\
	 data index 1 u 1:4 w linespoints pt 7 t "Err Prop_{M}"


reset 
set output "jack_Xi-C_vs_T.pdf"
data = "jack_Xi-C_vs_T.out"
set key top left
set grid
set logscale y
set xlabel "T"

Tc=2.28
set arrow from Tc, graph 0 to Tc, graph 1 nohead lc rgb "orange" lw 1.5

plot for [i=1:10] data index i u 1:2:4 w yerrorbars lc rgb "#e7298a" pt 7 notitle, \
     for [i=0:0] data index i u 1:2:4 w yerrorbars lc rgb "#e7298a" pt 7 t "{/Symbol C}", \
     for [i=1:10] data index i u 1:3:5 w yerrorbars lc rgb "#1b9e77" pt 5 notitle, \
     for [i=0:0] data index i u 1:3:5 w yerrorbars lc rgb "#1b9e77" pt 5 t "c"