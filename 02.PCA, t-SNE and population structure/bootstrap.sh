#!/bin/sh
export run=$1
export seed=$2
mkdir -p $run.Rep || exit 1
cp admixture.sh $run.Rep
cd $run.Rep
for i in {2..10}
do
bash admixture.sh $i $seed
done
