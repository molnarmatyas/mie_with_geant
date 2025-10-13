#!/bin/bash

#SET IT TO YOUR OWN!
HOME="/home/bporfy/Norma"

BUILD_DIR="$HOME/build"
MACRO_DIR="$HOME/macros"
INP_DIR="$HOME/poli_xsections"
NUM_TO_SHUF=$1
OUTPUT_NAME="$2"


# sanity check
echo "This is your home directory"
echo $HOME
sleep 5

# copy randomly selected files to the build directory
cp `ls -d "$INP_DIR/"* | shuf -n $NUM_TO_SHUF` "$BUILD_DIR/."


cd "$BUILD_DIR"

for file in limfo*.txt; do
  for i in $(seq 1 10);
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
GEANT_OUTPUT=$(ls -t ../build/output*.txt | head -n $NUM_TO_SHUF)

#run txtToHist on the outputs
j=0
for i in $GEANT_OUTPUT
do
  echo $i
  root -x -l -b -q txtToHist.C\(\"$i\",\"${OUTPUT_NAME}_${j}\"\)
  let j++
done
echo "Automated measurement done."

