#!/usr/bin/gnuplot -persist
filename = system("echo $FILENAME")

set size ratio -1
plot filename with linespoints linetype -1 linewidth 2 lc rgb "#6194B"