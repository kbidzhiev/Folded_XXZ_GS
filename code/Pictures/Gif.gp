set terminal gif animate delay 1
set output 'Sz_profile.gif'
stats 'Sz_profile.dat' nooutput
set xrange [-300:300]
set yrange [-0.6:0.6]
set ylabel "Sz"
set xlabel "x"
percentile="P5 P10 P20 P25 P50 P75"
cir=7;

f(x) = x;
do for [i=1:int(STATS_blocks)] {
    plot "Sz_profile.dat" index (i) u ($1):($2) w lp pt cir ps 1.5 lt rgb "red" title columnheader, "" index (i) u ((abs($1)<=f(i))? $1 : NaN ):(0.0) w l lw 10 ti "x=t"
}








