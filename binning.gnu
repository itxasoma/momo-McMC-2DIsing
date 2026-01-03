set term png
set output "binning_program_L_100_T_2.27.png
dades = "binning_program_L_100_T_2.27.dat"

set logscale x
set logscale y

plot dades u 1:3 t "E", dades u 1:5 t "M"

set xlabel 'Bin size m'
set ylabel 'Statistical error'

set key outside top center box horizontal

reset 
set output "jack.png"
set logscale x
dades= "timeseries_L100_T2.00_MCSTOT100000000 copy.dat"
set xlabel "MC steps"
set ylabel "Magnetization"
plot dades u 1:3 t "M", dades u 1:(abs($3)) t "|M|"