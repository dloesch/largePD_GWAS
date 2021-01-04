library(SeqArray)
library(SNPRelate)
library(gdsfmt)

source("/mnt/hdd/large-PD/v2/scripts/GENESIS/topmed.R")

argsv <- commandArgs(trailingOnly = TRUE)
#prefix 
prefix <- argsv[1]

#gdsfile
gdsfile <- argsv[2]


#seqarray
#gds <- seqOpen(gdsfile)

#snprelate
gds <- snpgdsOpen(gdsfile)

ibd <- snpgdsIBDKING(gds,type="KING-robust", autosome.only=FALSE, num.thread = countThreads())

outfile=paste(prefix,".ibd.RData",sep="")
save(ibd, file=outfile)

snpgdsClose(gds)
#seqClose(gds) #seqarray
