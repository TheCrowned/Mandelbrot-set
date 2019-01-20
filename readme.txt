mpic++ mandelbrot_par.cpp
mpirun -n 4 man_par 100 0.001 -1 -1 0.5 1 50 
sort -g -k2 -rk1 output.dat > sout
gnuplot 
plot "sout" with image
