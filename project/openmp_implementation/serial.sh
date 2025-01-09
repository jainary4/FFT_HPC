#!/bin/bash
EXECUTABLE="./fft"
OUTPUT_FILE="serial_results.txt"
# Clear the output file
echo "Serial Scaling Analysis Results" > $OUTPUT_FILE
echo -e "k\tN\tTime (s)" >> $OUTPUT_FILE

for k in {8..27}; do
  N=$((2**k))
  echo "Running for N = $N (k = $k)..."
   $EXECUTABLE $k 3 >> $OUTPUT_FILE

done
