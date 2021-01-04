library(SeqVarTools)
library(GENESIS)
library(gdsfmt)
library(SNPRelate) #choose if SeqArray is not optimized
#library(SeqArray)

source("/mnt/hdd/large-PD/v2/scripts/GENESIS/topmed.R")

argsv <- commandArgs(trailingOnly = TRUE)
gds_file <- argsv[1]
king_file <- argsv[2]
prefix <- argsv[3]
kin.type <- argsv[4]

#open gds file # if using seqarray format
#gds <- seqOpen(gds_file)
#if seqarray file
#seqData <- SeqVarData(gds)
#get sample.id
#sample.id <- seqGetData(gds, "sample.id")

#if using snprelate format
gds <- snpgdsOpen(gds_file)
sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))

type <- "IBD" #preferred method

if(type == "IBD"){
  load(king_file)
  sparse <- 0.01104854 #2^(-13/2), 5th degree
  king <- kingToMatrix(ibd, thresh=sparse)
}else{
  king <- king2mat(king_file, iids = sample.id)
}

save(king, file="kin.RData")

if (kin.type == "king") {
    kinMat <- king
    outfile <- paste(prefix,".pca.RData",sep="")
} else if (kin.type == "pcrelate") {
    load(paste(prefix,".pcrelate.RData",sep=""))
    kinMat <- pcrelateToMatrix(pcrel, scaleKin=2)
    outfile <- paste(prefix,".pca2.RData",sep="")
} else {
    stop("kinship method should be 'king' or 'pcrelate'")
}
message("Using ", kin.type, " kinship estimates")

kin_thresh <- 0.04419417 # 2^(-9/2), 3rd degree
n_pcs <- 20

#min maf and missing rate
maf=0.01
miss=0.01


#updated function
pca <- pcair(gds, kinobj=kinMat, divobj=king, num.cores=countThreads(),
             kin.thresh = kin_thresh, div.thresh= -kin_thresh, 
             maf=maf, missing.rate=miss, autosome.only=FALSE)

save(pca, file=outfile)

#seqClose(seqData)
snpgdsClose(gds)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
