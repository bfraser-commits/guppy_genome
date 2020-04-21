#!/bin/bash

# Load bedops
module load BEDOPS/2.4.26

# Coverage directory
cd $DIR

# Read the input and output from command line
WINDOW=$1

# Set pop_array as popmaps for each population males and females
pop_array=(GH_males GH_females GL_males GL_females LM_males LM_females UM_males UM_females LO_males LO_females UQ_males UQ_females)

# Loop over chrs & scafs
for CHR in $(cat ~/guppy_research/STAR/STAR.chromosomes.release.fasta.fai | cut -f1 | uniq)
do

echo "STARTING ${CHR}"

# Make reference
awk '{print $1,0,$2}' ~/guppy_research/STAR/STAR.chromosomes.release.fasta.fai | bedops --chop $WINDOW - | grep $CHR > winds_${CHR}_${WINDOW}.bed

# Loop over pops
for POP in "${pop_array[@]}"
do

# Set output
OUTPUT=${POP}_${CHR}_${WINDOW}_coverage

# Turn input into bed
awk -v CHR=$chr -v OFS="\t" 'NR>1 {print $1, $2, $2+50, 0, $3}' ${POP}_${CHR}_mashed2.bed > ${OUTPUT}_tmp.bed

# Do weighted mean within windows
bedmap --echo --wmean --skip-unmapped winds_${CHR}_${WINDOW}.bed ${OUTPUT}_tmp.bed | sed "s/|/\t/g" > ${OUTPUT}.bed

# Tidy
rm -f ${OUTPUT}_tmp.bed
done

# Tidy
rm -f winds_${CHR}_${WINDOW}.bed

done



