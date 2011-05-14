set terminal postscript  enhanced color eps
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set ytic auto                          # set ytics automatically
set title "Absolute error of Metropolis path tracing, when compared to Monte Carlo (after 1 mill. samples)"
set xlabel "X"
set ylabel "Error (W/m^2)"

set grid

plot "metropolis_error_1mill.dat" using (abs($2)) with lines lc rgb "red" title "Metropolis path tracing",\
"pathtracing_irrad.dat" using (abs($10-$10000)) with lines lc rgb "black" title "Monte carlo path tracing"
