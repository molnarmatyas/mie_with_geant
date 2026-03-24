#!/bin/bash

#SET IT TO YOUR OWN!
HOME="${PWD}" # should be where this script is located, and where the build and macro directories are located

BUILD_DIR="$HOME/build"
MACRO_DIR="$HOME/macros"
INP_DIR="/opt/cell_crosssections/mon/adda/high_fid"
NUM_TO_SHUF=$1
OUTPUT_NAME="$2"
RAND_POS_NUM=$3


# sanity check
echo "This is your home directory"
echo $HOME
sleep 5

# copy randomly selected files to the build directory
cp `ls -d "$INP_DIR/"* | shuf -n $NUM_TO_SHUF` "$BUILD_DIR/."


cd "$BUILD_DIR"
cmake ..
make -j2


#FIXME
for file in intensity_results_MON*.txt; do
  for i in $(seq 1 $RAND_POS_NUM);
  do
      x=$(seq -15 0.01 15 | shuf | head -n1)
      y=$(seq -15 0.01 15 | shuf | head -n1)
      echo $x,$y
      echo $PWD
      ./offset_runner.sh $x $x 1 $y $y 1 0 100000000 "$file"
  done
done


#select the last N outputs
cd "$MACRO_DIR"
GEANT_OUTPUT=$(ls -t ../build/output*.txt | head -n $(( $NUM_TO_SHUF*$RAND_POS_NUM )))

#run txtToHist on the outputs
j=0
for i in $GEANT_OUTPUT
do
  echo $i
  root -x -l -b -q txtToHist.C\(\"$i\",\"${OUTPUT_NAME}_${j}\"\)
  let j++
done
echo "Automated measurement done."

