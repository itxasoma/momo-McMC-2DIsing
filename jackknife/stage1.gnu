set terminal pdfcaire enhanced color

set output "block_stage1.pdf"
data = "jackerr_Xi-C_vs_block.out"

set xrange [2:]
set logscale xy
set grid
set key top left
set xlabel "block size"
set ylabel "{/Symbol s}"

plot data index 10 u 1:2 w linespoints pt 7 t "M", data index 0 u 1:3 w linespoints pt 7 t "E",\
	 data index 10 u 1:4 w linespoints pt 7 t "ErrProp(M)"