set term png 
set output "P3-20-21-ﬁg1.png"

set title "Gràfica 1"
set xlabel "E" 
set ylabel "D(E)"
set key outside 
set grid xtics
set grid ytics

plot [-1:7] [-200:400] "P3-20-21-res.dat" u 1:2 title "D(E)" w lp pointtype 7 pointsize 0.5 lt rgb "violet", "P3-20-21-res.dat" u 1:3 title "D'(E)" w lp pointtype 7 pointsize 0.5 lt rgb "blue"
