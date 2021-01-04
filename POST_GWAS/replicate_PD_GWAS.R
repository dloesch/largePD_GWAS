library(liftOver)
library(gwascat)
library(rtracklayer)



#replicate chang et al. 2017 and nalls 2014
nalls <- read.table("nalls_2014.GWAS_signif.txt", header=TRUE)
nalls$C <- paste0("chr",nalls$C)
chang <- read.delim("chang_etal.novel_loci.txt", header=TRUE)
chang$CHR <- paste0("chr", chang$CHR)

#liftover nalls data
path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
ch = import.chain(path)
temp <- makeGRangesFromDataFrame(nalls, keep.extra.columns = TRUE, start.field = "Position", end.field = "Position", seqnames.field = "C")
new <- liftOver(temp, chain=ch)
new = unlist(new)
genome(new) = "hg38"
new <- as.data.frame(new)
new <- new[-c(3:5)]
colnames(new)[1:2] <- c("CHR", "POS")
nalls <- new

#liftover chang data
temp <- makeGRangesFromDataFrame(chang, keep.extra.columns = TRUE, start.field = "BP", end.field = "BP", seqnames.field = "CHR")
new <- liftOver(temp, chain=ch)
new = unlist(new)
genome(new) = "hg38"
new <- as.data.frame(new)
new <- new[-c(3:5)]
colnames(new)[1:2] <- c("CHR", "POS")
chang <- new


#load in LARGE-PD GWAS results
assoc.file <- paste0("../GWAS/results/large.assoc_single.ALL.RData")
load(assoc.file)

#replicate nalls
snps <- nalls$SNP

rep_nalls <- nalls[c("SNP","Nearest_gene", "OR_Joint" )]

for(snp in snps){
  chr <- sub("chr","",nalls$CHR[nalls$SNP == snp])
  pos <- nalls$POS[nalls$SNP == snp]
  EA <- toupper(nalls$Effect_allele[nalls$SNP == snp])
  temp <- results[results$chr == chr & results$pos == pos,]
  if(nrow(temp) == 0){
    rep_nalls$BETA[rep_nalls$SNP == snp] <- NA
    rep_nalls$P[rep_nalls$SNP == snp] <- NA
  }else{
    if(temp$alt == EA){
      rep_nalls$BETA[rep_nalls$SNP == snp] <- temp$beta
    }else{
      rep_nalls$BETA[rep_nalls$SNP == snp] <- temp$beta*-1
    }
    rep_nalls$P[rep_nalls$SNP == snp] <- temp$Score.pval
  }
}

rep_nalls$OR <- exp(rep_nalls$BETA)

rep_nalls$CONCORDANCE <- ifelse((rep_nalls$OR_Joint >1 & rep_nalls$OR > 1) | (rep_nalls$OR_Joint < 1 & rep_nalls$OR < 1), TRUE, FALSE)
colnames(rep_nalls) <- c("SNP", "GENE", "OR_NALLS", "BETA_LARGE", "P_LARGE", "OR_LARGE", "CONCORDANCE")

#for chang
snps <- chang$SNP

rep_chang <- chang[c("SNP","Candidate_gene", "ORJoint" )]
rep_chang$ORJoint[17] <- chang$ORdiscovery[17]
for(snp in snps){
  chr <- sub("chr","",chang$CHR[chang$SNP == snp])
  pos <- chang$POS[chang$SNP == snp]
  EA <- toupper(chang$Effect_allele[chang$SNP == snp])
  temp <- results[results$chr == chr & results$pos == pos,]
  if(nrow(temp) == 0){
    rep_chang$BETA[rep_chang$SNP == snp] <- NA
    rep_chang$P[rep_chang$SNP == snp] <- NA
  }else{
    if(temp$alt == EA){
      rep_chang$BETA[rep_chang$SNP == snp] <- temp$beta
    }else{
      rep_chang$BETA[rep_chang$SNP == snp] <- temp$beta*-1
    }
    rep_chang$P[rep_chang$SNP == snp] <- temp$Score.pval
  }
}

rep_chang$OR <- exp(rep_chang$BETA)

rep_chang$CONCORDANCE <- ifelse((rep_chang$ORJoint >1 & rep_chang$OR > 1) | (rep_chang$ORJoint < 1 & rep_chang$OR < 1), TRUE, FALSE)

colnames(rep_chang) <- c("SNP", "GENE", "OR_CHANG", "BETA_LARGE", "P_LARGE", "OR_LARGE", "CONCORDANCE")

#foo et al. 

snps <- c("rs246814","rs9638616")

foo <- as.data.frame(snps)
foo$GENE <- c("SV3C", "WBSCR17")
foo$CHR <- c(5,7)
foo$POS <- c(76303383,71285507)
foo$Effect_allele <- c("T", "T") #the paper said rs246814's effect allele is C, but their given AF indicates it is actually the T
foo$OR_DISC <- c(1.11,1.14)
foo$OR_EUR <- c(1.07, 1.00)
foo$OR_ALL <- c(1.11, 1.02)

for(snp in snps){
  chr <- sub("chr","",foo$CHR[foo$snps == snp])
  pos <- foo$POS[foo$snps == snp]
  EA <- toupper(foo$Effect_allele[foo$snps == snp])
  temp <- results[results$chr == chr & results$pos == pos,]
  if(nrow(temp) == 0){
    foo$BETA[foo$snps == snp] <- NA
    foo$P[foo$snps == snp] <- NA
  }else{
    if(temp$alt == EA){
      foo$BETA[foo$snps == snp] <- temp$beta
    }else{
      foo$BETA[foo$snps == snp] <- temp$beta*-1
    }
    foo$P[foo$snps == snp] <- temp$Score.pval
  }
}

foo$OR <- exp(foo$BETA)
foo$CONCORDANCE <- ifelse((foo$OR_ALL >1 & foo$OR > 1) | (foo$OR_ALL < 1 & foo$OR < 1), TRUE, FALSE)

write.table(rep_chang, "large.chang_etal.replication.txt", sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)
write.table(rep_nalls, "large.nalls_etal_2014.replication.txt", quote=FALSE, row.names = FALSE, col.names = TRUE)
write.table(foo, "large.foo_etal.replication.txt", sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)
