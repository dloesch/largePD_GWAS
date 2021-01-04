library(SeqVarTools)
library(SNPRelate)
library(SeqArray)
source("/mnt/hdd/large-PD/v2/scripts/GENESIS/topmed.R")
source("/mnt/hdd/large-PD/v2/scripts/GENESIS/filterVariants.R")

#get prefix,gdsfile, and chr
argsv <- commandArgs(trailingOnly=T)
prefix <- argsv[1]
gdsfile <- argsv[2]

#snp characteristics
maf <- 0.01
miss <- 0.01

#create snprelate style gds file
new_gds <- paste0(prefix,".pruned.snp.gds")
seqGDS2SNP(gdsfile,new_gds)

gdsfile <- new_gds
gds <- snpgdsOpen(gdsfile)


outfile <- paste(prefix,".ALL.grm.RData", sep="")
  

grm <- snpgdsGRM(gds,method="GCTA", maf=maf, missing.rate= miss, 
                 autosome.only = FALSE, num.thread = countThreads())
  
save(grm,file=outfile)
  
snpgdsClose(gds)
