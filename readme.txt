mpic++ mandelbrot_par.cpp
mpirun -n 4 man_par 100 0.001 -1 -1 0.5 1 50 
sort -g -k2 -rk1 output.dat > sout
gnuplot 
plot "sout" with image

rows_per_task = 5 was around the best value, although it did not seem to yield significant speedup anyway. Tested with step 0.001, max_iterations 100, in the set {1 2 5 10 15 20 50 75 100 250 500 800} and run a couple times (file output disabled). Reasonable to believe the best depends on step size.
