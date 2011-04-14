set terminal png size 800, 600
set output "irrad_pathtracing.png"
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
#set boxwidth 0.9 absolute
set style fill solid 1.00 border -1
set style histogram clustered gap 1 title offset character 0, 0, 0
set datafile missing '-'
set style data histograms
set xtics border in scale 1,0.5 nomirror rotate by -70 offset character 0, 0, 0
set ytic auto                          # set ytics automatically
set title "Irradiance estimate using path tracing"
set xlabel "Samples"
set ylabel "Irradiance"
plot "irrad_pathtracing.dat" using 2:xtic(1) with lines title "Irradiance" lc rgb "black"
