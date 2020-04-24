##########################################################
# Fis calculation for chr12 for each population
# Scripts calculates Fis for 10kb windows
##########################################################

setwd("~/Documents/Genome_paper/")

lib<-c("pbapply","hierfstat","data.table","vcfR","adegenet","parallel","ggplot2","patchwork","cowplot")
lapply(lib,library,character.only=T)

# Run over all the pops
pops<-c("GH","GL","LM","UM","LO","UQ")
lapply(pops,function(pop){

# Read in data
vcf<-read.vcfR(paste0("~/Documents/Genome_paper/data/STAR_data/6_nations_STAR_hets_",pop,".vcf.gz"))

# Pops
all_females<-read.table("~/Documents/Genome_paper/data/six_nations_females.popmap")[,1]
inds<-colnames(vcf@gt)[2:ncol(vcf@gt)]
females<-inds[inds %in% all_females]
males<-inds[!(inds %in% all_females)]
popmap<-rep("male",length(inds))
popmap[inds %in% females]<-"female"

# Winds
windows<-seq(0,max(as.integer(vcf@fix[,2])),10000)
windows2<-windows+10000


 wind_Fis<-mclapply(1:length(windows),function(x){
   print(x)
   
   if(length(which(as.integer(vcf@fix[,2]) > windows[x] & 
                   as.integer(vcf@fix[,2]) < windows2[x])) < 2){
     out<-data.frame(chr="chr12",
                     start=windows[x],
                     end=windows2[x],
                     Fis=NA)
     return(out)  
   } else {
   
   # Subset VCF
   vcf_sub<-vcf[which(as.integer(vcf@fix[,2]) > windows[x] & 
                  as.integer(vcf@fix[,2]) < windows2[x]),]
   
   # Convert
   dat<-vcfR2genind(vcf_sub)
   
   # Make custom popmap
   popmap2<-popmap[inds %in% rownames(dat$tab)]
   pop(dat)<-popmap2
   dat2<-genind2hierfstat(dat)
   
   # Calculate
   stats<-basic.stats(dat2,diploid = 2,digits = 2)
   
   # Out
   out<-data.frame(chr="chr12",
                   start=windows[x],
                   end=windows2[x],
                   Fis=mean(na.omit(stats$perloc$Fis)))
   return(out)
   }
 #})))
},mc.cores=detectCores()-1)

saveRDS(wind_Fis,file = paste0("outputs/",pop,"_chr12_fis_calcs.rds"))
#wind_Fis<-readRDS(paste0("outputs/",pop,"_chr12_fis_calcs.rds"))
fis_dd<-data.frame(rbindlist(wind_Fis))

write.table(fis_dd,
			paste0("outputs/",pop,"_chr12_hets_Fis_windows.txt"),
			row.names=F,quote=F,sep="\t")

})
