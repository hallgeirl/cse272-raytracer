set terminal postscript enhanced color eps
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set ytic auto                          # set ytics automatically
set title "Mean square error of Metropolis path tracing with a sphere at (0,0,0), compared to Monte Carlo"
set xlabel "Number of sample rays"
set ylabel "Mean square error"
set logscale y

set grid

plot "metropolis_sphere_msq.dat" with lines lc rgb "red" title "Metropolis path tracing, mean square error"
