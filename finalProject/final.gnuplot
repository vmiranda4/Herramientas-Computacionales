set terminal postscript color enhanced
set output "strange_attractor.eps"

set xlabel "eje x" 
set ylabel "eje y"
set zlabel "eje z"
set title "Grafica strange attractor"
set grid
splot "coordinates.txt" using 2:3:4 with lines title "xyz" ,\