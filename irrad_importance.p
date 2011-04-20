set terminal png size 1600, 1200
set output "irrad_importance.png"
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set ytic auto                          # set ytics automatically
set title "Irradiance estimate using importance sampling"
set xlabel "Samples"
set ylabel "Irradiance"
plot [0:1000000] "irrad_importance.dat" every 100 with lines title "Irradiance" lc rgb "black"
