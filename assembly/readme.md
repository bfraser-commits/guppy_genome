# GitHub Pull

```
git clone https://github.com/bfraser-commits/guppy_genome.git code
````

# Falcon Assembly


```
ln -s code/falcon-assembly.cfg
ln -s code/reads-fq.fofn
fc_run falcon-assembly.cfg
```

# Falcon Unzip


```
ln -s code/reads-bam.fofn
ln -s code/falcon-unzip.cfg
fc_unzip.py falcon-unzip.cfg

```

# Quiver Polishing

```
mkdir 4-polish
cd 4-polish
cat ../3-unzip/*.fa > contigs.fasta
samtools faidx contigs.fasta
pbalign reads-bam.fofn contigs.fasta contigs.bam
quiver contigs.bam -r contigs.fasta -o contigs_quiver.fasta --parametersSpec P6-C4
grep ">" contigs_quiver.fasta | grep "_" |cut -f2 -d ">" | SeqFilter --ids - contigs_quiver.fasta --out all_h_ctg.fasta
grep ">" contigs_quiver.fasta | grep "_" |cut -f2 -d ">" | SeqFilter --ids - contigs_quiver.fasta --ids-exclude --out all_p_ctg.fasta
cd ..
```

# Falcon Phase

```
mkdir 5-phase
cd 5-phase
ln -s ../code/falcon-phase.json contig.json
ln -s ../4-polish/all*.fasta ./
sed 's/\|quiver//' all_p_ctg.fasta > all_p_ctg.clean.fasta
sed 's/\|quiver//' all_h_ctg.fasta > all_h_ctg.clean.fasta  
scrub_names.pl all_p_ctg.fasta  all_h_ctg.fasta > name_mapping.txt
snakemake -s snakefile --verbose -p
cd..
```

# Pilon Polishing

```
mkdir 6-polish
cd 6-polish
ln -s ../5-phase/PRET.phase_0.fasta
bwa index PRET.phase_0.fasta
bwa mem PRET.phase_0.fasta $PE1 $PE2 | samtools view -bS - | samtools sort > PRET.phase_0.bam
samtools index PRET.phase_0.bam
java -Xmx200G -jar pilon-1.16.jar --genome PRET.phase_0.fasta --frags PRET.phase_0.bam --output PRET.phase_0.polished.fasta
cd..
```

# Genetic Mapping

```
mkdir 7-linkage
cd 7-linkage
cat ../6-polished/PRET.phase_0.polished.fasta | sed '/>/s/\|pilon//g' > genome.fasta
ln -s ../code/weights_*.txt ./
bwa index genome.fasta
bwa mem genome.fasta GM1.fasta | samtools view -bS -q 60 | samtools sort | bedtools bam2bed -i - > GM1.bed
bwa mem genome.fasta GM2.fasta | samtools view -bS -q 60 | samtools sort | bedtools bam2bed -i - > GM2.bed
ln -s weights_split.txt weights.txt
python -m jcvi.assembly.allmaps merge GM1.bed GM2.bed -o map.bed
python -m jcvi.assembly.allmaps split map.bed --chunk=4 > breakpoints.bed
# Additional breaks from Hi-C analysis where added at this point. See breakpoint.bed in repo for final edits!
sed -e '/^>/! s/[[:lower:]]/N/g' genome.fasta > genome_n.fasta
python -m jcvi.formats.fasta gaps genome_n.fasta
python -m jcvi.formats.sizes agp genome.fasta
python -m jcvi.formats.agp mask genome.fasta.agp breakpoints.genome.refined.bed --splitobject --splitsingle
python -m jcvi.formats.agp build genome.fasta.masked.agp genome.fasta genome.SPLIT.fasta
bwa index genome.SPLIT.fasta
bwa mem genome.SPLIT.fasta GM1.fasta | samtools view -bS -q 60 | samtools sort | bedtools bam2bed -i - > GM1.SPLIT.bed
bwa mem genome.SPLIT.fasta GM2.fasta | samtools view -bS -q 60 | samtools sort | bedtools bam2bed -i - > GM2.SPLIT.bed
bwa mem genome.SPLIT.fasta GM3.fasta | samtools view -bS -q 60 | samtools sort | bedtools bam2bed -i - > GM3.SPLIT.bed
bwa mem genome.SPLIT.fasta GM4.fasta | samtools view -bS -q 60 | samtools sort | bedtools bam2bed -i - > GM4.SPLIT.bed
rm weights.txt
ln -s weights_map.txt weights.txt
python -m jcvi.assembly.allmaps merge GM1.SPLIT.bed GM2.SPLIT.bed GM3.SPLIT.bed GM4.SPLIT.bed -o release.bed
python -m jcvi.assembly.allmaps path release.bed genome.SPLIT.fasta
cd ..
ln -s 7_linkage/release.fasta
ln -s 7_linkage/release.agp
```
