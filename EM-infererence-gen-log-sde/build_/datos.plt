 set terminal wxt size 971,600 enhanced font 'Verdana,10' persist
 set title "function g"
 set xlabel "m"
 set ylabel "g(m)"
 set key outside 
 set arrow from 20, graph 0 to 20, graph 1 nohead
 plot "datos.dat" using 1:2 with lines title "X_t:[1]"
