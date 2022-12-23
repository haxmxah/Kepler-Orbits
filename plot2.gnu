set term png 
set output "P3-20-21-ﬁg2.png"

set title "Gràfica 2"
set xlabel "E" 
set ylabel "D(E)"
set key outside 
set grid xtics
set grid ytics

plot "P3-20-21-res.dat" i 2 u 2:3 title "x(E)" w lp pointtype 7 pointsize 0.5 lt rgb "violet"
