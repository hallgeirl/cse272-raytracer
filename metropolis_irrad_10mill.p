set terminal postscript  enhanced color eps
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set ytic auto                          # set ytics automatically
set title "Irradiance, Metropolis path tracing at 10 mill. samples"
set xlabel "X"
set ylabel "Irradiance (W/m^2)"

set grid

plot "metropolis_irrad.dat" using 10000 with lines lc rgb "red" title "Metropolis path tracing",\
"pathtracing_irrad.dat" using 100 with lines lc rgb "black" title "Monte carlo path tracing"
