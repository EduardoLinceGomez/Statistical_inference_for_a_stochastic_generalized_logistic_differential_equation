# define line styles using explicit rgbcolor names
#
set style line 1  linetype 1 linecolor rgb "black"  linewidth 1.000 \
    pointtype 6 pointsize default pointinterval 0
set style line 2  linetype 2 linecolor rgb "orange"  linewidth 2.000 \
    pointtype 2 
# setting terminal
#set term wxt size 971,600
set terminal pngcairo size 971,600 enhanced font 'Verdana,9'
set output 'plot_g_m.png'
# svg
#set terminal svg size 971,600 fname 'Verdana, Helvetica, Arial, sans-serif' \
#fsize '9' rounded dashed
#set output 'plot_g_m.svg'
set title "function g"
set xlabel "m"
set ylabel "g(m)"
x_position = 2

set arrow 1 at x_position, graph 0 to x_position, graph 1 nohead lc "red" dt 4
set arrow 2 from 0,0 to 5,0 nohead
set key outside
#set arrow from 0.20, graph 0 to 0.20, graph 1 nohead
plot "datos.dat" using 1:2 ls 1 title "X_t:[1]"


