cd build
#time ./exampleB1 run2.mac -t 1 &> example0.log
time ./exampleB1 run2.mac -t 2> mie.err 1> mie.out
cd ..
grep "mag,phi,theta:" build/mie.out | awk '{print $(NF-1), $NF}' > build/mie_filter.out
#root.exe -b -q 'make_histogram.C("build/output100M.root", 100, "#Theta", "pngname")'
root.exe -b -q 'make_histogram_old.C("build/mie_filter.out", 100, "#theta", "thetadist.png")'
