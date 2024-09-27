#!/usr/bin/gnuplot -persist
filename = system("echo $FILENAME")

#set xrange[-1.0:28.0]
#set yrange[-1.0:1.0]
set size ratio -1
plot filename with linespoints linetype -1 linewidth 2 lc rgb "#6194B"