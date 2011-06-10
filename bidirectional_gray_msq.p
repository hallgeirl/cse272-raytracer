set terminal postscript enhanced color eps
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set ytic auto                          # set ytics automatically
set title "Mean square error of bidirectional Metropolis path tracing, compared to Monte Carlo"
set xlabel "Number of sample rays"
set ylabel "Mean square error"
set logscale y

set grid

plot "bidirectional_gray_msq.dat" with lines lc rgb "red" title "Bidirectional Metropolis path tracing, mean square error"
