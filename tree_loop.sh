#!/bin/bash

for i in {1..10}
do
  
  echo "Submitting script $i"
  qsub ~/path/to/tree.sh $i 64 25 &

done

exit 0