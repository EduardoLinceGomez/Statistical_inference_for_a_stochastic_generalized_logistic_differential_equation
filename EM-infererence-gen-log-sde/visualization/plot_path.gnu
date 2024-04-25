 # define line styles using explicit rgbcolor names
#
set style line 1  linetype 1 linecolor rgb "black"  linewidth 1.000 \
    pointtype 6 pointsize default pointinterval 0
set style line 2  linetype 2 linecolor rgb "orange"  linewidth 2.000 \
    pointtype 2
set style line 3  linetype 3 linecolor rgb "black"  linewidth 1.000 \
    pointtype 0 
# setting terminal
# setting terminal
#set term wxt size 971,600
set terminal pngcairo size 971,600 enhanced font 'Verdana,9'
set output 'plot_X_t.png'
# sv
#set terminal svg size 971,600 fname 'Verdana, Helvetica, Arial, sans-serif' \
#fsize '9' rounded dashed
#set output 'plot_g_m.svg'
#set terminal wxt size 971,600 enhanced font 'Verdana,10' persist
 set title "Sampled paths"
 set xlabel "t"
 set ylabel "X_t"
 set key outside 
 plot "tray.dat" using 1:2 ls 3 title "X_t[1]" 
