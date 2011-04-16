set terminal png size 1600, 1200
set output "irrad_pathtracing.png"
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set ytic auto                          # set ytics automatically
set title "Irradiance estimate using path tracing"
set xlabel "Samples"
set ylabel "Irradiance"
plot [10000:10000000] "irrad_pathtracing.dat" every 500 with lines title "Irradiance" lc rgb "black"
