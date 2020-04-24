##################################################################
# Script for running PopGen data for Guppy Genome Paper - Fraser et al. 2020
# Calculates male/female Fst, dxy and pi for each population
# Calculates these for 1kb and 10kb windows
##################################################################

# Set
setwd("/Users/jw962/Documents/Genome_paper")

# Get libs
lib<-as.vector(c("pbapply","PopGenome","ggplotify","cowplot","gridExtra","viridis","grid","data.table","HiTC","Sushi","ggplot2","parallel","dplyr","tidyr"))
lapply(lib,library,character.only=TRUE)

vcf<-"data/STAR_data/6_nations_STAR_GATK4_pop_SNP.gatk.bi.miss.maf.final.filtered.recode.vcf.gz"

# Catch errors
wind_sizes=c(1000,10000)

for(i in 1:length(wind_sizes)){

# Set window size
	wind_size<-wind_sizes[i]


# Get length
chr_length <- read.table("data/genomes/STAR/STAR.chromosomes.release.fasta.fai", header=F)
chrs<-as.character(read.table("data/STAR_data/STAR_chrs.txt")[,1])

# Run over chrs
all_chr<-data.frame(rbindlist(pblapply(1:length(chrs),function(chr){
#all_chr<-data.frame(rbindlist(lapply(141:length(chrs),function(chr){
    
# chr name  
tid = chrs[chr]
print(tid)
# chr length
topos = chr_length[chr_length$V1 == tid,]$V2


if(wind_size < topos )
{
  
# Read in the VCF
snp <- readVCF(file = vcf, tid=tid, frompos = 1, topos = topos, numcols=1000000, include.unknown=TRUE)

# Set the blokes and sheilas
GH_M = c("GH11","GH12","GH13","GH14","GH15","GH17","GH18","GH19","GH20")
GH_F = c("GH1","GH10","GH2","GH3","GH4","GH5","GH6","GH7","GH8","GH9")
  
MH_M = c("LM1","LM2","LM4","LM5","LM3","LM10","LM6","LM7","LM8","LM9")
MH_F = c("LM11","LM12","LM13","LM14","LM15","LM16","LM17", "LM18","LM19","LM20")
 
OH_M <- c("LO2","LO7","LO10","LO5","LO3","LO4","LO6","LO12","LO11","LO12")
OH_F <- c("LO1" ,"LO13","LO15","LO16","LO17","LO18","LO19","LO20","LO21","LO22")
 
GL_M = c("GL11", "GL12", "GL13", "GL14", "GL15", "GL16", "GL17", "GL18", "GL19", "GL20")
GL_F =c("GL1", "GL10", "GL2", "GL5", "GL6", "GL7", "GL8", "GL9")
 
ML_M = c("UM2","UM3","UM6","UM1","UM5","UM7")
ML_F = c("UM10","UM11","UM12","UM13","UM14","UM15","UM16","UM17","UM18","UM8","UM9")
 
OL_M <- c("UQ1","UQ5","UQ6","UQ7","UQ8","UQ2","UQ3","UQ4")
OL_F <- c("UQ10","UQ11","UQ12","UQ13", "UQ14", "UQ15","UQ16","UQ17","UQ18", "UQ19", "UQ20","UQ9")

GH_list<-list(GH_M,GH_F)
GL_list<-list(GL_M,GL_F)
MH_list<-list(MH_M,MH_F)
ML_list<-list(ML_M,ML_F)
OH_list<-list(OH_M,OH_F)
OL_list<-list(OL_M,OL_F)

comparison_lists<-list(GH_list,GL_list,MH_list,ML_list,OH_list,OL_list)
names<-c("GH","GL","MH","ML","OH","OL")
 
# Lapply over list of comparisons
chr_popgen<-data.frame(rbindlist(mclapply(1:length(comparison_lists),function(x){

snp<-set.populations(snp,do.call("list",comparison_lists[[x]]),diploid=TRUE)
snp<-set.outgroup(snp,FALSE)

# Split into windows
win_SNP_10k<-sliding.window.transform(snp,width=wind_size,jump=wind_size,type=2)

#do pop stats
win_SNP_10k <-F_ST.stats(win_SNP_10k,mode="nucleotide")

win_SNP_10k <-neutrality.stats(win_SNP_10k,FAST=FALSE)

win_SNP_10k <-diversity.stats.between(win_SNP_10k,nucleotide.mode=TRUE)

# Get centre of the window
 genome.pos_10K <- sapply(win_SNP_10k@region.names, function(x){
    split <- strsplit(x," ")[[1]][c(1,3)]
   val   <- mean(as.numeric(split))
    return(val)
  })
   
# Output results matrix
PG_out<-data.frame(chr = tid,
                   FST=win_SNP_10k@nucleotide.F_ST,
                   n.seg.M=win_SNP_10k@n.segregating.sites[,1],
                   n.seg.F=win_SNP_10k@n.segregating.sites[,2],
                   dxy=(win_SNP_10k@nuc.diversity.between/wind_size), 
                   pi.M=(win_SNP_10k@nuc.diversity.within/wind_size)[,1],
                   pi.F=(win_SNP_10k@nuc.diversity.within/wind_size)[,2],
                   start=genome.pos_10K -(wind_size/2)-1.5,
                   end=genome.pos_10K +(wind_size/2)-1.5,
                   pop=names[x])
colnames(PG_out)[c(2,5)]<-c("FST","DXY")
#PG_out<-na.omit(PG_out)
PG_out[PG_out$FST < 0,"FST"]<-0
rownames(PG_out)<-NULL

return(PG_out)
},mc.cores=3)))

write.table(na.omit(chr_popgen),
            paste0("outputs/popgen_outputs/6nations_STAR_popgen_",wind_size,"_",tid,".txt"),row.names=F,quote=F,sep="\t")
return(chr_popgen)
}
})))
}




