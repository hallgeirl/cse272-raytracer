set terminal postscript  enhanced color eps
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set ytic auto                          # set ytics automatically
set title "Irradiance, Adaptive photon mapping at 10 mill. samples"
set xlabel "X"
set ylabel "Irradiance (W/m^2)"

set yrange [0:0.8]
set grid

plot "adaptiveppm_irrad.dat" using 100 with lines lc rgb "red" title "Adaptive photon mapping",\
"pathtracing_irrad.dat" using 100 with lines lc rgb "black" title "Monte Carlo path tracing"
