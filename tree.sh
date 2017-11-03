#!/bin/bash
#proteus header

#load r open

#whatever

WORK_DIR=~/Documents/rcrust

i=$1
nc=$2
seed=$3

R CMD BATCH "--args $i $nc $seed" $WORK_DIR/proteus_tree.R $WORK_DIR/proteus_tree.out &

exit 0