set term png 
set output "P3-20-21-ﬁg3.png"

set title "Gràfica 3"
set xlabel "E" 
set ylabel "dF(E),df_{approx}(E)"
set key outside 
set grid xtics
set grid ytics


plot [0:2.5] [-15:10] "P3-20-21-res3-n34.dat" u 1:3 title "dF(E)" w l,"P3-20-21-res3-n34.dat" u 1:4 title "dF_{approx}(E)"
