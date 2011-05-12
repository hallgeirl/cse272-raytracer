set terminal postscript  enhanced color eps
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set ytic auto                          # set ytics automatically
set title "Irradiance estimate using path tracing"
set xlabel "Samples (x1000)"
set ylabel "X"
set zlabel "Irradiance"

set ytics ("-1" 0, "-0.5" 25, "0" 50, "0.5" 75, "1" 100)

#set pm3d
#set dgrid3d 100,100,2
set grid
set view 45,15

splot "pathtracing_irrad.dat" matrix every 10:1 with lines lw 0.1 lc rgb "black" title ""
