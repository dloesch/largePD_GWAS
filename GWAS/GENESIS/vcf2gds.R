library(SeqArray)
library(tools)
library(SNPRelate)


source("/mnt/hdd/large-PD/v2/scripts/GENESIS/topmed.R")

#get vcf and gds file
argsv <- commandArgs(trailingOnly=T)

vcffile <- argsv[1]
gdsfile <- argsv[2]

#bcf file
isBCF <- file_ext(vcffile) == "bcf"
if (isBCF) {
  ## use bcftools to read text
  vcffile <- pipe(paste("bcftools view", vcffile), "rt")
}
## pick format fields to import
fmt.import <- "DS"

#set reference genome
ref <- "hg38"

## write to the scratch disk of each node
#gdsfile.tmp <- tempfile()
#message("gds temporarily located at ", gdsfile.tmp)

###this function does not work on 3.5.1 on EDN cluster
seqVCF2GDS(vcffile, gdsfile,reference=ref, 
           storage.option="LZMA_RA", parallel=TRUE, scenario = "imputation")
#set parallel=countThreads() if able to multithread

##alternative function using snprelate format
#snpgdsVCF2GDS(vcffile,gdsfile, ignore.chr.prefix = "chr")

#convert to seqarray format, if desired
#seqSNP2GDS(gdsfile.tmp, gdsfile)
#outfile <- paste("SEQ.",gdsfile, sep="")
#seqSNP2GDS(gdsfile, outfile)

## copy it
#file.copy(gdsfile.tmp, gdsfile)
## remove the tmp file
#file.remove(gdsfile.tmp)
