~
All work was done on a laptop with Ubuntu 18.04; Intel Core i3-7100U @ 2.40 GHz x 4, 8GB of DDR4 RAM and SSD of 256 GB. 
(Some, in particular plotted ones) Measures were taken on Cheetah cluster. However, due to *my mistake* in submitting the job (not using a nodelist), those measures only used 4 cores and did not exploit parallelism well.

~
Usage instructions are provided in the compiled versions of programs. 
Anyway, a sample invocation would be:

man_seq <lower left X coord> <lower left Y coord> <upper right X coord> <upper right Y coord> <max iterations> <step size> (<rows per task>)

where:
- the first 4 parameters define the grid in which computation is run. Always tested with (-1, -1) and (0.5, 1), as that's where the core of the set is
- <max iterations> is the max number of iterations after which a point is considered to be divergent (i.e. not belonging to mandelbrot set)
- <step size> controls the granularity of computation
- <rows per task>, only present in the parallel version, controls how many rows each task is made of. This is useful in particular for large workloads (i.e. small step sizes), for which sending just one line to be processed has a considerable overhead vs sending 5-10 of them.

~
To obtain a plot of the mandelbrot set, make sure output lines are not commented out in the source and use something like the following:

mpic++ mandelbrot_par.cpp
mpirun -n 4 man_par -1 -1 0.5 1 100 0.001 50 
sort -g -k2 -rk1 output.dat > sout
gnuplot 
plot "sout" with image

An example can be found as mandelbrot-set.png. The sequential/parallel versions do not differ in generating the image, so the sequential can be used as well. Keep in mind that a step size smaller than 0.001 (I tried with 0.0001) does not work because the output file is several GBs in size, and gnuplot goes out of memory while trying to plot it.

~
During my tests, I found that rows_per_task = 5 was around the best value, although it did not seem to yield a *significant* speedup anyway. Tested with step 0.001, max_iterations 100, in the set {1 2 5 10 15 20 50 75 100 250 500 800} and run a couple times (with file output disabled). Measures results can be found in measure/times-rows_per_task.txt.
I find reasonable to believe that the best value depends on the step size: the smaller the step, the more a bigger rows_per_task could help lowering the communication overhead and exploiting parallelism.

~
As expected, I found that max_iterations does not make much difference whether the set is computed in parallel or sequential version, as they are a serial task (no parallelization happens on discovering *when/whether* a points diverges).

~
I obtained my plot of the speedup curve of sequential vs parallel with

set logscale xy
plot "measure/times-merge.txt" using 2:4 with linespoints ls 1 title "Parallel", "measure/times.txt" using 2:5 with linespoints ls 2 title "Sequential"

this can be found as speedup.png.

The speedup, albeit scarce, shows that both programs scale in exponential complexity, but the parallel version is quite faster than the sequential one.

~
Debug mode can be enabled in both sequential and parallel versions by tweaking the value of the constant DEBUG at the top.
