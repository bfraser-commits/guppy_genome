#!/bin/bash

# Load what we need, including the deeptools anaconda library
module load Anaconda3/5.2.0
source activate deeptools-env

# Make output dir
mkdir $MASTER/outputs/coverage
mkdir $MASTER/outputs/coverage/$DATASET
mkdir $MASTER/outputs/coverage/$DATASET/$POP
OUT=$MASTER/outputs/coverage/$DATASET/$POP

# Make the 'blacklisted' regions of 100kb at the beginning and end of chr
CHR_END=$(grep -w "chr${MOAB_JOBARRAYINDEX}" $MASTER/data/STAR.chromosomes.release.fasta.fai | cut -f2)
CHR_END2=$(bc <<< "scale = 5; $CHR_END - 100000")
echo -e "chr${MOAB_JOBARRAYINDEX}\t0\t100000\tNA\nchr${MOAB_JOBARRAYINDEX}\t$CHR_END2\t$CHR_END\tNA" > $OUT/chr${MOAB_JOBARRAYINDEX}_black.bed

# Make a list of scaffolds to exclude for normalisation calculations
scafs=(`cat ~/guppy_research/STAR/STAR.chromosomes.release.fasta.fai | cut -f1 | grep -v "chr"`)
scaf_blacks=""i
for i in "${scafs[@]}"
do
    scaf_blacks+="$i "
done

# Loop over bams
# Note - The Effective Genome Size has been reduced by 100,000 to account for Blacklist
for bam in $(ls ${BAM_DIR}/*.bam | grep $POP)
do
#bamCoverage -b $bam -bl $OUT/chr${MOAB_JOBARRAYINDEX}_black.bed -o ${bam}.chr${MOAB_JOBARRAYINDEX}.coverage.bed -p 4 --outFileFormat=bedgraph --binSize=100 --smoothLength=300 --effectiveGenomeSize=528351853 -ignore $scaf_blacks --region=chr${MOAB_JOBARRAYINDEX} --normalizeUsing=RPGC
bamCoverage -b $bam -bl $OUT/chr${MOAB_JOBARRAYINDEX}_black.bed -o ${bam}.chr${MOAB_JOBARRAYINDEX}.coverage.bed -p 4 --outFileFormat=bedgraph --binSize=50 --smoothLength=75 --effectiveGenomeSize=528351853 --region=chr${MOAB_JOBARRAYINDEX} --normalizeUsing=RPGC
#bamCoverage -b $bam -bl $OUT/chr${MOAB_JOBARRAYINDEX}_black.bed -o ${bam}.chr${MOAB_JOBARRAYINDEX}.coverage.bed -p 4 --outFileFormat=bedgraph --binSize=50 --smoothLength=75 --effectiveGenomeSize=707997498 --region=chr${MOAB_JOBARRAYINDEX} --normalizeUsing=RPGC
mv ${bam}.chr${MOAB_JOBARRAYINDEX}.coverage.bed $OUT
done

# Tidy up
rm -f $OUT/chr${MOAB_JOBARRAYINDEX}_black.bed

# Exit environment
source deactivate
