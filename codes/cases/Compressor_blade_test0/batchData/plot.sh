#!/usr/bin/gnuplot -persist

plot '0.dat' with linespoints linetype -1 linewidth 2 lc rgb "#6194B",'1.dat' with linespoints linetype -1 linewidth 2 lc rgb "#f58231",'2.dat' with linespoints linetype -1 linewidth 2 lc rgb "#ffe119",'3.dat' with linespoints linetype -1 linewidth 2 lc rgb "#bfef45",'4.dat' with linespoints linetype -1 linewidth 2 lc rgb "#3cb44b",'5.dat' with linespoints linetype -1 linewidth 2 lc rgb "#42d4f4",'6.dat' with linespoints linetype -1 linewidth 2 lc rgb "#4363d8",'7.dat' with linespoints linetype -1 linewidth 2 lc rgb "#911eb4",'8.dat' with linespoints linetype -1 linewidth 2 lc rgb "#f032e6",'9.dat' with linespoints linetype -1 linewidth 2 lc rgb "#000075",'10.dat' with linespoints linetype -1 linewidth 2 lc rgb "#aaffc3",

while(1){
	replot
	pause 1
}
