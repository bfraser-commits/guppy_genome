###########################################################################################################################################################################
# This script takes a directory of deeptools coverage bams merged over populations and averaged into windows and examines coverage and coverage ratios between HP and LP pops
###########################################################################################################################################################################

## Load packages
lib = c("parallel","data.table","tidyr","dplyr","ggplot2")
lapply(lib, library, character.only=TRUE)

# Get parameters from commandline to find outputs
setwd("~/deeptools")
input<-"~/deeptools/outputs/coverage/six_nations_STAR_FINAL/pop_avgs_filtered"
  
# Define pops
pops<-c("GH_males","GH_females",
"GL_males","GL_females",
"LM_males","LM_females",
"UM_males","UM_females",
"LO_males","LO_females",
"UQ_males","UQ_females")

pairs<-lapply(seq(1,length(pops),2),function(x){return(pops[c(x,x+1)])})
names(pairs)<-c("GH","GL","MH","ML","OH","OL")
rivers<-rep(c("Guanapo","Marianne","Oropouche"),each=4)
pred<-rep(c("High","High","Low","Low"),3)

# Define window size
window_vec<-c("1000","10000")

mclapply(window_vec,function(wind){

WINDOW=wind

# List all the files
inputs<-list.files(input)

# Read in all the data to 1 file
dd<-data.frame(rbindlist(lapply(1:length(pops),function(x){

# Read in chr/scf at a time
chr_scfs<-as.character(read.table("data/STAR.chromosomes.release.fasta.fai")[,1])
chr_scfs<-chr_scfs[grep("alt",chr_scfs,invert=T)]

# Get file
to_read<-grep(pops[x],inputs,value=T)
to_read<-grep(paste0("_",WINDOW,"_"),to_read,value=T)

tmp_dd<-data.frame(rbindlist(lapply(1:length(to_read),function(y){

# Read file
tmp<-data.frame(fread(paste0(input,"/",to_read[y])))

# Rename
colnames(tmp)<-c("chr","BP1","BP2","coverage")

# Add cols
tmp$pop<-pops[x]

if(x %in% seq(1,length(pops),2)){
tmp$sex<-"Male"
} else {
tmp$sex<-"Female"
}

tmp$river<-rivers[x]
tmp$pred<-pred[x]

return(tmp)
})))

return(tmp_dd)
})))

######################################################################
# Do Plots and visualisations
chr_scfs<-as.character(read.table("data/STAR.chromosomes.release.fasta.fai")[,1])
chr_scfs<-chr_scfs[grep("alt",chr_scfs,invert=T)]

plotting_function<-lapply(chr_scfs[1:23],function(x){
print(x)

tmp<-dd[dd$chr == x,]

if(nrow(tmp) > 0){

# First do plot of raw data
raw_plot<-ggplot(tmp,aes(x=BP1,y=coverage,colour=sex,linetype=sex))+
          geom_line()+
          facet_grid(river~pred,scales="free_y")+
          xlab(x)+
          scale_x_continuous(breaks = seq(0,signif(max(tmp$BP2),2),by=2000000),
                       labels=seq(0,signif(max(tmp$BP2),2),by=2000000)/1000000)+
                          theme_bw()+
          theme(panel.grid = element_blank(),
          axis.title = element_text(size=24),
          axis.text = element_text(size=14),
          strip.text = element_text(size=12),
          legend.position="right")+
          ylab("Coverage (log)")+
          labs(colour="Sex",linetype="Sex")+
          scale_colour_manual(values=c("red2","blue2"))


# Now we need to calculate coverage ratio per pair
ratios<-data.frame(rbindlist(lapply(1:length(pairs),function(y){

tmp2<-tmp[tmp$pop %in% pairs[[y]],]
BP_to_keep<-duplicated(tmp2$BP1)
tmp2<-tmp2[tmp2$BP1 %in% tmp2$BP1[BP_to_keep],]

M_cov<-tmp2[tmp2$pop == pairs[[y]][1],"coverage"]
F_cov<-tmp2[tmp2$pop == pairs[[y]][2],"coverage"]

                
                out<-data.frame(chr=rep(tmp2$chr),
                BP1=tmp2[tmp2$pop== pairs[[y]][1],"BP1"],
                BP2=tmp2[tmp2$pop== pairs[[y]][1],"BP2"],
                pop=rep(names(pairs)[y]),
                river=tmp2[tmp2$pop== pairs[[y]][1],"river"],
                pred=tmp2[tmp2$pop== pairs[[y]][1],"pred"],
                cov_ratio=M_cov/F_cov,
                M_cov_est=M_cov,
                F_cov_est=F_cov)
                
# Correct for cov_ratio "Infs"
out$cov_ratio_corrected<-out$cov_ratio
out[out$cov_ratio_corrected == "Inf","cov_ratio_corrected"]<-max(out[out$cov_ratio != "Inf","cov_ratio"])

# Remove NAs (which come from 0s for each)
out<-na.omit(out)

# Correct for 0 values and set as minimum
out[out$cov_ratio_corrected == 0,"cov_ratio_corrected"]<-min(out[out$cov_ratio != 0,"cov_ratio"])
                
# Also estimate a spline
spline_mean <- smooth.spline(y = log2(out$cov_ratio_corrected), x=out$BP1, spar = 0.05)
out$ratio_spline<-spline_mean$y

spline_mean2 <- smooth.spline(y = out$M_cov_est, x=out$BP1, spar = 0.05)
out$M_spline<-spline_mean2$y

spline_mean3 <- smooth.spline(y = out$F_cov_est, x=out$BP1, spar = 0.05)
out$F_spline<-spline_mean3$y
   
return(out)
})))

# Plot raw-estimate splines
          spline_raw_plot<-ggplot(ratios,aes(x=BP1,y=M_spline))+
          geom_line(colour="blue2",linetype="dashed")+
          geom_line(data=ratios,aes(x=BP1,y=F_spline),colour="red2")+
          facet_grid(river~pred,scales="free_y")+
          xlab(x)+
          scale_x_continuous(breaks = seq(0,signif(max(ratios$BP2),2),by=2000000),
                       labels=seq(0,signif(max(ratios$BP2),2),by=2000000)/1000000)+
                          theme_bw()+
          theme(panel.grid = element_blank(),
          axis.title = element_text(size=24),
          axis.text = element_text(size=14),
          strip.text = element_text(size=12),
          legend.position="none")+
          ylab("Raw Coverage Splines")

# Plot raw values as above as log2-transform
raw_ratio_plot<-ggplot(ratios,aes(x=BP1,y=log2(cov_ratio)))+
          geom_line()+
          facet_grid(river~pred,scales="free_y")+
          xlab(x)+
          scale_x_continuous(breaks = seq(0,signif(max(ratios$BP2),2),by=2000000),
                       labels=seq(0,signif(max(ratios$BP2),2),by=2000000)/1000000)+
                          theme_bw()+
          theme(panel.grid = element_blank(),
          axis.title = element_text(size=24),
          axis.text = element_text(size=14),
          strip.text = element_text(size=12),
          legend.position="none")+
          ylab("Coverage Ratio Male:Female")
          
          spline_ratio_plot<-ggplot(ratios,aes(x=BP1,y=ratio_spline))+
          geom_line()+
          facet_grid(river~pred,scales="free_y")+
          xlab(x)+
          scale_x_continuous(breaks = seq(0,signif(max(ratios$BP2),2),by=2000000),
                       labels=seq(0,signif(max(ratios$BP2),2),by=2000000)/1000000)+
                          theme_bw()+
          theme(panel.grid = element_blank(),
          axis.title = element_text(size=24),
          axis.text = element_text(size=14),
          strip.text = element_text(size=12),
          legend.position="none")+
          ylab("Coverage Ratio Male:Female")

# Also plot a smoothed spline

# Plot them both together
pdf(paste0("figs/STAR_",x,"_Sex_",WINDOW,"_Coverage_Raw_figs.pdf"),width=16)
print(raw_plot)
print(spline_raw_plot)
print(raw_ratio_plot)
print(spline_ratio_plot)
dev.off()

return(ratios)
}
})

# Now just get the actual table for outputting
final_out<-lapply(chr_scfs,function(x){

tmp<-dd[dd$chr == x,]
print(x)

if(nrow(tmp) > 0){

# Now we need to calculate coverage ratio per pair
ratios<-data.frame(rbindlist(lapply(1:length(pairs),function(y){

tmp2<-tmp[tmp$pop %in% pairs[[y]],]
BP_to_keep<-duplicated(tmp2$BP1)
tmp2<-tmp2[tmp2$BP1 %in% tmp2$BP1[BP_to_keep],]

if(nrow(tmp2) > 0){

M_cov<-tmp2[tmp2$pop == pairs[[y]][1],"coverage"]
F_cov<-tmp2[tmp2$pop == pairs[[y]][2],"coverage"]

out<-data.frame(chr=rep(tmp2$chr),
                BP1=tmp2[tmp2$pop== pairs[[y]][1],"BP1"],
                BP2=tmp2[tmp2$pop== pairs[[y]][1],"BP2"],
                pop=rep(names(pairs)[y]),
                river=tmp2[tmp2$pop== pairs[[y]][1],"river"],
                pred=tmp2[tmp2$pop== pairs[[y]][1],"pred"],
                cov_ratio=M_cov/F_cov,
                M_cov_est=M_cov,
                F_cov_est=F_cov)
                
return(out)
}
})))
return(ratios)
}
})

# Save the ratios
write.table(data.frame(rbindlist(final_out)),
            paste0("outputs/STAR_coverage_ratios_",WINDOW,".txt"),row.names=F,quote=F,sep="\t")
            
},mc.cores=3)
