#!/bin/bash

# Make output directory
mkdir $MASTER/$DATASET/$POP/pop_avg_filtered

cd $MASTER/$DATASET/$POP/

# This version of the script includes a filtering step to remove all rows with coverage > 4
parallel -j 12 "cat *.chr{}.* | awk '$4 < 4 {print $0}' | datamash -H -s -g 2 mean 4 | sort -nk1 > ${POP}_chr{}_mashed.bed" ::: {1..23}

# Loop over all chr and scafs
for scaf in $(cat /gpfs/ts0/home/jw962/HP_LP_trials/deeptools/data/STAR.chromosomes.release.fasta.fai | cut -f1 | uniq)
do
cat *.${scaf}.* | awk '$4 < 4 {print $0}' | datamash -H -s -g 2 mean 4 | sort -nk1 > ${POP}_${scaf}_mashed.bed
awk -v chr="$scaf" '{print chr,"\t",$0}' ${POP}_${scaf}_mashed.bed > ${POP}_${scaf}_mashed2.bed
done
rm -f *_mashed.bed
mv *mashed* pop_avg_filtered
done
