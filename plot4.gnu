set term png 
set output "P3-20-21-ﬁg4.png"

set title "Gràfica 4"
set xlabel "E" 
set ylabel "F(E)"
set key outside 
set grid xtics
set grid ytics




plot [0:2.5] [-15:10] "P3-20-21-res3-n420.dat" u 1:3 title "dF(E)" w l,"P3-20-21-res3-n420.dat" u 1:4 title "dF_{approx}(E)"
