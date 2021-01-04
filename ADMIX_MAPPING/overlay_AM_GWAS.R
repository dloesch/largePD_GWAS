library(data.table)
library(coda)
library(liftOver)
library(gwascat)
library(rtracklayer)
library(readxl)


#load joint test
load("./Admixture_mapping/Admixture_mapping_LARGE_PD_EUR.NAM_031220.RData")
chrs <- unique(EUR.NAM$chr[EUR.NAM$Joint.pval < 5E-04]) #just run on chromosomes with "interesting" results

#set significance level using empirical autoregression
tests <- c()

for(chr in 1:22){
  foo <- EUR.NAM[EUR.NAM$chr == chr,]
  fit <- spectrum0.ar(foo$Joint.pval)
  tests <- c(tests, fit$spec)
}
total <- sum(tests)
signif <- 0.05/total
print(signif)

#take whichever is more conservative
signif <- ifelse(signif > 5E-05, 5E-05, signif)


##prepare dataframe that will have both AM summary stats and finemapping stats
finemap <- as.data.frame(chrs)
colnames(finemap) <- "CHR"

#pdf("large-PD.AM_finemap.pdf")

for(chr in chrs){
  
  #the admixure mapping results from Andrea did not have the physical position, so I used 1000 Genomes data to get position information
  
  #get hg19 positions and rs ids from 1000 Genomes file
  kgfile <- paste0("./1000_genomes/plink/EUR/EUR.chr",chr, ".bim")
  kg <- fread(kgfile, header=FALSE, stringsAsFactors = FALSE, data.table=FALSE)
  kg <- kg[c("V1","V2", "V4")]
  colnames(kg) <- c("CHR", "SNP", "POS")
  kg$CHR <- paste0("chr", kg$CHR)
  
  
  #lift KG positions to hg38 to match GWAS data
  path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
  ch = import.chain(path)
  temp <- makeGRangesFromDataFrame(kg, keep.extra.columns = TRUE, start.field = "POS", end.field = "POS")
  new <- liftOver(temp, chain=ch)
  new = unlist(new)
  genome(new) = "hg38"
  new <- as.data.frame(new)
  new <- new[c("SNP", "start")]
  colnames(new)[2] <- "POS"
  
  
  #merge KG with AM results
  foo <- EUR.NAM[EUR.NAM$chr == chr,]
  foo <- merge(foo, new, by.x="snpID", by.y="SNP", all.x=TRUE)
  foo <- foo[complete.cases(foo),]
  
  
  #ID boundaries of the "main" peak
  pvals <- unique(foo$Joint.pval[foo$chr == chr])
  pvals <- pvals[order(pvals)][1:2]
  start <- min(foo$POS[foo$chr == chr & foo$Joint.pval %in% pvals])
  end <-  max(foo$POS[foo$chr == chr & foo$Joint.pval %in% pvals])
  
  #subset AM results to area +/- 1E06 of AM peak for plotting
  foo <- foo[foo$POS > start - 1E06 & foo$POS < end + 1E06,]
  smooth_locus = smooth.spline(x=foo$POS,y=-log10(foo$Joint.pval), spar=0.1)
  
  #get single variant results from LARGE-PD file and find intersection with AM peak
  results.file <- paste0("/mnt/hdd/large-PD/v2/GWAS/results/large.assoc_single.chr",chr,".RData")
  load(results.file)
  assoc <- assoc[assoc$pos > start - 1E06 & assoc$pos < end + 1E06,]
  
  
  plot(-log10(assoc$Score.pval)~assoc$pos,type="p", pch=16, xlab="position (bp)", ylab="-log10 p-value", 
       main=paste0("Chromosome ",chr," locus"), cex=0.9) #ylim = c(0, xmax)
  lines(smooth_locus, col="dodgerblue", lwd=3)
  
  
  #####FINEMAP PROCEDURE#######
  
  #save AM summary statistics
  finemap$PEAK[finemap$CHR == chr] <- paste0(start, "-", end)
  finemap$P[finemap$CHR == chr] <- min(foo$Joint.pval)
  finemap$P_ADJ[finemap$CHR == chr] <- min(foo$Joint.pval)*(0.05/signif)
  
  #subset association data
  assoc <- assoc[assoc$pos > start & assoc$pos < end,]
  assoc <- assoc[assoc$freq > 0.01 & assoc$freq < 0.99,]
  
  #add top variant in peak
  top <- assoc[assoc$Score.pval == min(assoc$Score.pval),]
  finemap$TOP_SNP[finemap$CHR == chr] <- new$SNP[new$POS == top$pos]
  finemap$TOP_SNP_POS[finemap$CHR == chr] <- top$pos
  finemap$TOP_SNP_PVAL[finemap$CHR == chr] <- top$Score.pval
  finemap$TOP_SNP_P_ADJ[finemap$CHR == chr]  <- top$Score.pval*nrow(assoc)
  
  
  ##highlight gene in AM peak###
  
  if(chr == 6){
    gene <- c(166409364,166906451)
    temp <- assoc[assoc$pos > gene[1] & assoc$pos < gene[2],]
    finemap$GENE[finemap$CHR == chr] <- "RPS6KA2"
    
    #add points to figure
    points(temp$pos, -log10(temp$Score.pval),col="firebrick2", pch=16, cex=0.9)
    par(font=4)
    legend(median(temp$pos),1, legend= "RPS6KA2", xjust=0.5, cex=0.8, adj =0.15, bg="white")
  }
  
  
  if(chr == 14){
    stxbp6 <- c(24809654,25050297)
    temp <- assoc[assoc$pos > stxbp6[1] & assoc$pos < stxbp6[2],]
    finemap$GENE[finemap$CHR == chr] <- "STXBP6"
    
    #add points to figure
    points(temp$pos, -log10(temp$Score.pval),col="red", pch=16, cex=0.9)
    par(font=4)
    legend(median(temp$pos),1, legend= "STXBP6", xjust=0.5, cex=0.8, adj =0.15, bg="white")
  }
  
  
  if(chr == 17){
    finemap$GENE[finemap$CHR == chr] <- NA
    
  }
  
  if(chr == 21){
    gene <- c(44885949,44931989)
    temp <- assoc[assoc$pos > gene[1] & assoc$pos < gene[2],]
    finemap$GENE[finemap$CHR == chr] <- "ITGB2"
    
    #add points to figure
    points(temp$pos, -log10(temp$Score.pval),col="firebrick2", pch=16, cex=0.9)
    par(font=4)
    legend(median(temp$pos),1, legend= "ITGB2", xjust=0.5, cex=0.8, adj =0.25, bg="white")
  }
  
  #add single ancestry results to dataframe
  load("/mnt/hdd/large-PD/v2/AM/Admixture_mapping/Admixture_mapping_LARGE_PD_AFR_031220.RData")
  load("/mnt/hdd/large-PD/v2/AM/Admixture_mapping/Admixture_mapping_LARGE_PD_NAM_031220.RData")
  
  AFR <- merge(AFR, new, by.x="snpID", by.y="SNP", all.x=FALSE, all.y=FALSE)
  AFR <- AFR[AFR$POS > start - 1E06 & AFR$POS < end + 1E06,]
  best_afr <- min(AFR$pval)
  
  NAM <- merge(NAM, new, by.x="snpID", by.y="SNP", all.x=FALSE, all.y=FALSE)
  NAM <- NAM[NAM$POS > start - 1E06 & NAM$POS < end + 1E06,]
  best_nam <- min(NAM$pval)
  
  
  #save single ancestry results
  finemap$AFR_PVAL[finemap$CHR == chr] <- best_afr
  finemap$ADJ_AFR_PVAL[finemap$CHR == chr] <- best_afr*(0.05/signif)
  finemap$NAM_PVAL[finemap$CHR == chr] <- best_nam
  finemap$ADJ_NAM_PVAL[finemap$CHR == chr] <- best_nam*(0.05/signif)
}

#dev.off()


#set adjusted pvalues to 1 if > 1
finemap$P_ADJ <- ifelse(finemap$P_ADJ > 1, 1, finemap$P_ADJ)
finemap$TOP_SNP_P_ADJ <- ifelse(finemap$TOP_SNP_P_ADJ > 1, 1, finemap$TOP_SNP_P_ADJ)
finemap$ADJ_AFR_PVAL <- ifelse(finemap$ADJ_AFR_PVAL > 1, 1, finemap$ADJ_AFR_PVAL)
finemap$ADJ_NAM_PVAL <- ifelse(finemap$ADJ_NAM_PVAL > 1, 1, finemap$ADJ_NAM_PVAL)

write.table(finemap, "large-PD.AM.finemap.updated.txt", sep="\t", row.names = FALSE, quote=FALSE, col.names = TRUE)
