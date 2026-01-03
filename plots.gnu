set term png
set out 'E_L_100_T_2_MC_106.png'

dades = 'timeseries_L100_T2.00_MCSTOT106.dat'
set xlabel 'MC time'
set ylabel 'Energy'

set logscale x
unset key
plot dades u 1:2 w lp 

set out 'M_L_100_T_2_MC_106.png'
set ylabel 'Magnetization'
plot dades u 1:3 w lp
