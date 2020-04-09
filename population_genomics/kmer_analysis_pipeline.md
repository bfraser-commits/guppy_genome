# Kmer analysis pipeline 


### Concatenate clean read files 
```
GL_females=(GL1 GL2 GL5 GL6 GL7)
GL_males=(GL11 GL12 GL13 GL14 GL15)
GH_females=(GH1 GH2 GH3 GH4 GH5)
GH_males=(GH11 GH12 GH13 GH14 GH15)
ML_females=(UM8 UM9 UM10 UM11 UM12)
ML_males=(UM1 UM2 UM6 UM7 UM5)
MH_females=(LM11 LM12 LM13 LM14 LM15)
MH_males=(LM1 LM2 LM3 LM4 LM5)
OL_females=(UQ9 UQ10 UQ11 UQ12 UQ13)
OL_males=(UQ1 UQ5 UQ6 UQ7 UQ8)
OH_females=(LO15 LO17 LO18 LO19 LO20)
OH_males=(LO2 LO3 LO4 LO5 LO6)
```

```
for i in ${GL_females[@]}; do zcat "$i_1.fq.gz" > GL_females_1.fq.gz; done
for i in ${GL_females[@]}; do zcat "$i_2.fq.gz" > GL_females_2.fq.gz; done
```

etc ...

### Run lighter for read correction for each populations
`lighter -r left.fq -r right.fq -k 31 730000000`

### Run jellyfish for kmerising reads

Count kmers

`jellyfish count -m 31 -s 10000000 -t 8 -C input`

Merge kmer binaries

`jellyfish merge -o output.jf output\_*`

Dump

`jellyfish dump -c sample.jf | sort -k 1 > counts.tsv`

####  Now we have a counts.tsv file representing the k-mer abundances for males and females from each pop ... 


### Use DE-kupl joinCounts to join all k-mer abundance files into a single matrix

`joinCounts GL_females.tsv GH_females.tsv ML_females.tsv MH_females.tsv OL_females.tsv OH_females.tsv \
GL_males.tsv GH_males.tsv ML_males.tsv MH_males.tsv OL_males.tsv OH_males.tsv > raw-counts-matrix.tsv`

### Create female filtered matrix:
<= 1 in females and >=1 males

_NB format of the matrix file looks like this_

`kmer_ID	GL_F	GH_F	ML_F	MH_F	OL_F	OH_F	GL_M	GH_M	ML_M	MH_M	OL_M	OH_M`

Calculate mean and upper 95% of the male columns using awk and then use these values in python script:

## Run kmer_picking_read1.py to pick read1 kmers

### When kmer picking has finished you may need to remove duplicate reads (kmers likely to appear in same read multiple times)

`sed -e '/^>/s/$/@/' -e 's/^>/#/' filein.fa | tr -d '\n' | tr "#" "\n" | tr "@" "\t" | sort -u -k1,1 | sed -e 's/^/>/' -e 's/\t/\n/' > fileout.fa`

Then remove first > at top of file
... and grep out just the header names:

`cat read1_ymers.fasta | grep "^>"  > $PATH/kmer_picking_outputs/tsv_files/ymers.tsv `

_use this tsv file to pick read 2s using the pick read 2 script)_

## Run kmer_associate_read2.py to get associated read2s for the read1 kmers

When you have the read 2's, you have to reorder them so they are in the same order as the read1s:

`sed -e '/^>/s/$/@/' -e 's/^>/#/' filein.fa | tr -d '\n' | tr "#" "\n" | tr "@" "\t" | sort -u -k1,1 | sed -e 's/^/>/' -e 's/\t/\n/' > fileout.fa`

and again remove the > from the first line

### Assemble read1 and read2 using Spades:

`spades.py -1 read1_ymers.fasta -2 read2_ymers.fasta \
--only-assembler --careful -o $dir -t 8`

### Align to genome:
`minimap2 --cs -ax sr $genome --secondary=no ymers_scaffolds.fasta > output.ymer.bam`




