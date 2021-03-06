# _Poecilia reticulata_ male genome annotation pipeline


## What we need

Genome file:

`genome=$PATH/guppy_genome/chromosome.release.fasta`

RNAseq data:

`rna=$PATH/rnaseq_data`

Sharma data:
		
`sharma=$rna/sharma`
		
Ovary and testes data:

`gonad=$rna/gonad`
		
 Working directory:

`working=$PATH/annotation`

##  Step 1: Genome masking

`cd $working/masking`

`ln -s $genome .`

### Export path to red and run
```
export PATH=/data/programs/RED/redUnix64:/data/programs/RED/redmask-master:$PATH
ingenome=chromosome.release.fasta
outgenome=chromosome.release
redmask.py -i $ingenome -o $outgenome
```

Masked genome: chromosome.release.softmasked.fa
num scaffolds : 267
assembly size : 746,481,348 bp
masked repeats : 210,899,717 bp (28.25%)

Set masked genome as new genome variable for mapping

`genome=$working/masking/chromosome.release.softmasked.fa`

##  Step 2a: Sharma paired-end read mapping

`cd $working/sharma`

Get paths to clean fastq reads and symlink them in:

`find $sharma -name "*.fq" > sharma.fq.rnaseq.files`

`while read p; do ln -s $p .; done < sharma.fq.rnaseq.files`

Get file headers
`ls *fq | cut -f1 -d '_'| uniq > PE_headers`

### Map using hisat2

`export PATH=/data/programs/hisat2-master/:$PATH`

`hisat2-build $genome ${genome%.fasta}_Hi2index`

### run hisat2 in parallel
dry run:

`cat PE_headers | parallel --dryrun 'hisat2 -p 8 -x ${genome%.fasta}_Hi2index -S {}.sam -q -1 {}rna_R.1clean.fq -2 {}rna_R.2clean.fq '`

run:

`cat PE_headers | parallel 'hisat2 -p 8 -x ${genome%.fasta}_Hi2index -S {}.sam -q -1 {}_1.fastq.ctq.fq -2 {}_2.fastq.ctq.fq'`

convert to BAM:

`parallel "samtools view -bS -@ 4 {.}.sam > {.}.bam" ::: *.sam`

sort and index:

`parallel "samtools sort {} {.}.sorted.bam; samtools index {.}.sorted.bam" ::: *.bam`

sanity check on the bam files:

`samtools quickcheck -v *.sorted.bam > bad_bams.fofn   && echo 'all ok' || echo 'some files failed check, see bad_bams.fofn'`

Clean up:
```
find . -type l -exec unlink {} \;
find . -type f -name "*.sam" rm {} \;
find . -type f ! -name '*.sorted.bam' -delete
```

##  Step 2b: Gonad single-end read mapping

`cd $working/gonad`

Get paths to clean fastq reads and symlink them in:

```
find $gonad -name "*.fq" > gonad.fq.rnaseq.files
while read p; do ln -s $p .; done < gonad.fq.rnaseq.files
```

get RNAseq file headers:

`ls *fq | cut -f1 -d '_'| uniq > single_RNAseq`

PATH to hisat2:

`export PATH=/data/programs/hisat2-master/:$PATH`

### run hisat2 mapping:

`cat single_RNAseq | parallel 'hisat2 -p 8 -x ${genome%.fasta}_Hi2index -S {.}.sam -q -U {}'`

convert to BAM:

`parallel "samtools view -bS -@ 4 {.}.sam > {.}.bam" ::: *.sam`
sort and index:

`parallel "samtools sort -m 1G -@ 5 {} -o {.}.sorted.bam ; samtools index {.}.sorted.bam" ::: *.bam`
sanity check:

`samtools quickcheck -v *sorted.bam > bad_bams.fofn   && echo 'all ok' || echo 'some files failed check, see bad_bams.fofn'`

Clean up:
```
find . -type l -exec unlink {} \;
find . -type f -name "*.sam" rm {} \;
find . -type f ! -name '*.sorted.bam' -delete
```

##  Step 3: BRAKER2 annotation 

### All the dependency paths
```
export PATH=/data/programs/augustus-3.3.1/bin/:/data/programs/BRAKER.v.2.1.1/scripts/:$PATH
export GENEMARK_PATH=/data/programs/gm_et_linux_64/gmes_petap
export AUGUSTUS_CONFIG_PATH=/data/programs/Augustus-3.3.1/config
export AUGUSTUS_SCRIPTS_PATH=/data/programs/Augustus-3.3.1/scripts
export BAMTOOLS_PATH=/data/programs/bamtools-2.5.1/bin/
export BLAST_PATH=/data/programs/ncbi-blast-2.2.31+/bin
export SAMTOOLS_PATH=/data/programs/samtools-1.7/
export ALIGNMENT_TOOL_PATH=/data/programs/gth-1.7.0-Linux_x86_64-64bit/bin
```

`cd $working/braker`

get a list of sorted bam files

`find $working -name "*sorted.bam" > bam.files`

make a soft link to all of them

```
while read p; do ln -s $p .; done < bam.files
ls -1 *sorted.bam | sed 's/$/,/' | tr "\n" " " | sed 's/ //g' > sorted_bam_file_names
bamfiles=sorted_bam_file_names
```

### Run BRAKER

`braker.pl --species=guppy_male --genome=$genome --bam=$bamfiles --cores 48`

##  Step 4: Quantitative assessment of the annotation

the path for the braker output is here:

`braker=$working/braker`

assess GFF

`gff_file=/$braker/augustus.hints.gtf`

reference

`reference=$working/masking/chromosome.release.softmasked.fa`

### Get CDS and protein fastas

`gffread "$gff_file" -g "$reference" -x "$reference"_cds.fa -y "$reference"_prot.fa`
  
working directory

`OHR=$working/OHR_assessment`

### Use xiphophorus_maculatus as model species proteome

```
mkdir xiphophorus_maculatus
cd xiphophorus_maculatus
wget ftp://ftp.ensembl.org/pub/release-96/fasta/xiphophorus_maculatus/pep/Xiphophorus_maculatus.X_maculatus-5.0-male.pep.all.fa.gz
gunzip Xiphophorus_maculatus.X_maculatus-5.0-male.pep.all.fa.gz
```

### cd hit clustering

cd-hit at 90% similarity

```
infile=Xiphophorus_maculatus.X_maculatus-5.0-male.pep.all.fa
outfile=Xiphophorus_maculatus.X_maculatus-5.0-male.pep.all.90.fa
cd-hit -i $infile -o $outfile -c 0.9 -T 16 -aS 0.8 -M 0
```

`ref_protein_set=Xiphophorus_maculatus.X_maculatus-5.0-male.pep.all.90.fa`
`query_protein_set=chromosome.release.softmasked.fa_prot.fa`

#### make blast db of genome set

`makeblastdb -dbtype prot -in $ref_protein_set`

### blast using parallel

`cat $query_protein_set | parallel --block 10k -k --recstart '>' --pipe "blastp -evalue 0.00001 -outfmt '6 qseqid sseqid qlen length slen qstart qend sstart send evalue bitscore pident' -db $ref_protein_set -query - " > $query_protein_set.v.$ref_protein_set.tsv`

`file=chromosome.release.softmasked.fa_prot.fa.v.Xiphophorus_maculatus.X_maculatus-5.0-male.pep.all.90.fa.tsv`
    
### filter on OHR > 0.95

`awk -F"\t" '{if ($6 > 0.95 ) print $0}' $file | wc -l` 
18265

get IDs of OHR_95

`awk -F"\t" '{if ($6 > 0.95 ) print $0}' $file | cut -f1 | sed '1d' > OHR_longest_gt_95`

```
reference=chromosome.release.softmasked.fa_prot.fa
fasta_header=OHR_longest_gt_95
samtools faidx $fasta_file
```

search and extract

`while read p; do samtools faidx $fasta_file $p ; done < $fasta_header > $fasta_header.fa`

### Then load results into eggnog and download

At an e-value of 1e-10 you get 15,263 genes


