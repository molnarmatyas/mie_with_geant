#!/bin/bash

BUILD_DIR="build"

# Copy all txt files to the build directory
cp /home/molnarmatyas/normawork/scattnlay-master/tests/python/limfo/limfo*.txt "$BUILD_DIR/"


cd "$BUILD_DIR"

for file in limfo*.txt; do
  for i in $(seq 1 20);
  do
      x=$(seq -5 0.01 5 | shuf | head -n1)
      y=$(seq -5 0.01 5 | shuf | head -n1)
      echo $x,$y
      ./offset_runner.sh $x $x 1 $y $y 1 0 100000000 "$file"
  done
done

