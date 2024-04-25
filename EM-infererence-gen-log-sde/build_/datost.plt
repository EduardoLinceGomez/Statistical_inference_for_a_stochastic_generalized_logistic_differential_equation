 set terminal wxt size 971,600 enhanced font 'Verdana,10' persist
 set title "Sampled paths"
 set xlabel "t"
 set ylabel "X_t"
 set key outside 
 plot "tray.dat" using 1:2 with lines title "X_t[1]" 
