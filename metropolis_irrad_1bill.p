set terminal postscript  enhanced color eps
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set ytic auto                          # set ytics automatically
set title "Irradiance estimate using metropolis sampling, 1 billion rays"
set xlabel "X"
set ylabel "Irradiance"
set yrange [0:0.8]

#set pm3d
#set dgrid3d 100,100,2
set grid

plot "metropolis_irrad.dat" using 1000000 with lines lc rgb "red" title "Metropolis path tracing", \
"pathtracing_irrad.dat" using 10000 with lines lc rgb "black" title "Monte Carlo path tracing (1 billion rays)"
