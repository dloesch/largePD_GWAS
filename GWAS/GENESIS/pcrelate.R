library(SeqVarTools)
library(GENESIS)
library(gdsfmt)
library(GWASTools)
source("/mnt/hdd/large-PD/v2/scripts/GENESIS/topmed.R")

argsv <- commandArgs(trailingOnly = TRUE)
#get prefix
prefix <- argsv[1]

#get pca file
pca_file <- paste(prefix,".pca.RData", sep="")
pca <- getobj(pca_file)


#get gdsfile
gdsfile <- argsv[2]

#if using snpgds format
#gds <- GdsGenotypeReader(gdsfile)
#genoData <- GenotypeData(gds)

#if using seqarray format
gds <- seqOpen(gdsfile)
seqData <- SeqVarData(gds)

variant.id <- NULL

var_block_size <- 1024
sample_block_size <- 10000
n_pcs <- 3
n_pcs <- min(n_pcs, length(pca$unrels))

#get sample.id
#sample.id <- seqGetData(gds, "sample.id") #seqarray
#sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))

#Min minor allele frequency
maf=0.01 #default of pcrelate function


#old pcrelate function
#pcrelate(seqData,
#         pcMat=pca$vectors[,1:n_pcs],
#         training.set=pca$unrels,
#         scan.include=sample.id, snp.include=variant.id,
#         write.to.gds=TRUE, gds.prefix=prefix,
#         scan.block.size=sample_block_size)


#iterator if using seqarray format
iterator <- SeqVarBlockIterator(seqData, variantBlock=sample_block_size)

#iterator if using snprelate format
#iterator <- GenotypeBlockIterator(genoData, snpBlock=var_block_size)


#updated pcrelate function
pcrel <- pcrelate(iterator,
                  pcs=pca$vectors[,1:n_pcs],
                  training.set=pca$unrels,
                  maf.thresh=maf)

outfile <- paste(prefix,".pcrelate.RData", sep="")
save(pcrel, file=outfile)

#close(iterator)
#for seqarray:
seqClose(seqData)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
