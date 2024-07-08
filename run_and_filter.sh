#!/bin/bash

cd build
./exampleB1 run2.mac 2> mie.err 1> mie.out
cd ..
grep "mag,phi,theta:" build/mie.out | awk '{print $(NF-1), $NF}' > mie.out
