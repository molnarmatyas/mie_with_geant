cd build
time ./mieNorma run2.mac -t 1 &> mieNorma.log
#time ./mieNorma run2.mac -t 2> mie.err 1> mie.out
cd ..
#grep "mag,phi,theta:" build/mie.out | awk '{print $(NF-1), $NF}' > build/mie_filter.out
root.exe -b -q 'plotCrossSection.C("build/output0.root", 100, "#theta", "thetadist.png",0.,3.14)' &> build/histprintout.log
#root.exe -b -q 'plotCrossSection.C("build/output100M.root", 100, "#theta", "thetadist.png",0.,3.14)'
#root.exe -b -q 'plotCrossSection_old.C("build/mie_filter.out", 100, "#theta", "thetadist.png",0.,3.14)'
