set terminal postscript  enhanced color eps
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set ytic auto                          # set ytics automatically
set title "Irradiance, Progressive photon mapping at 100 mill. samples"
set xlabel "X"
set ylabel "Irradiance (W/m^2)"

set yrange [0:0.8]

set grid

plot "progressive_irrad.dat" using 1000 with lines lc rgb "red" title ""
