library(SeqVarTools)
library(Biobase)
library(GENESIS)
library(SeqArray)
#library(GWASTools)
#library(SNPRelate)

source("/mnt/hdd/large-PD/v2/scripts/GENESIS/topmed.R")

#get chr
argsv <- commandArgs(trailingOnly=T)
chr <- argsv[1]

#get gds file
gdsfile <- argsv[2]

#get prefix
prefix <- argsv[3]

#get file names
null_model_file <- paste(prefix,".null_model.RData", sep="")
phenotype_file <- paste(prefix,".pheno.RData", sep="")


#analysis characteristics
out_prefix = paste(prefix,".assoc_single", sep="")
test_type= "score"
variant_block_size= 1024
pass_only=FALSE #keep as false if using seqarray format
maf_threshold=0.01
mac_threshold=NA #MAC filtering not supported with snpgds format


## open gds file

#seqArray
gds <- seqOpen(gdsfile)

#snprelate format
#gds <- snpgdsOpen(gdsfile)


# createSeqVarData object
# get phenotypes
load(phenotype_file)


annot <- phen$annot
pData(annot)$sample.id <- as.character(pData(annot)$sample.id) #LRAGE-PD has numeric ids, change to character
#if genotype data has more subjects than phenotype data, use the next 3 lines
s <- as.data.frame(seqGetData(gds,"sample.id"),stringsAsFactors=FALSE)
colnames(s) <- "sample.id"
pData(annot) <- merge(s, pData(annot),by.x="sample.id",by.y="sample.id",all.x=TRUE,all.y=TRUE, sort=FALSE)
pData(annot) <- merge(s, pData(annot),by.x="sample.id",by.y="sample.id",all.x=TRUE,all.y=TRUE, sort=FALSE)
seqData <- SeqVarData(gds, sampleData=annot)


# get null model
nullModel <- getobj(null_model_file)

# get samples included in null model
sample.id <- nullModel$sample.id


#ony works for seqarray format
if (as.logical(pass_only)) {
    filterByPass(seqData)
}

## MAC/MAF filtering

mac.min <- as.numeric(mac_threshold)
maf.min <- as.numeric(maf_threshold)

#filtering for snprelate file
#snpset <- snpgdsSelectSNP(gds, maf=maf.min, missing.rate=0.05, autosome.only = F)
#filtered <- unlist(unname(snpset))

#close snprelate gds file
#snpgdsClose(gds)

#the following code works for seqarray format gds files
if (!is.na(mac.min)) {
    filterByMAC(seqData, sample.id, mac.min=mac.min)
} 
#
if (!is.na(maf.min)) {
   filterByMAF(seqData, sample.id, maf.min=maf.min)
}

checkSelectedVariants(seqData)

# create iterator, seqarray
block.size <- as.integer(variant_block_size)
iterator <- SeqVarBlockIterator(seqData, variantBlock=block.size)


#iterator from GWASTools
#wants ID column to be called "scanID" for some reason
#pheno <- pData(phen$annot)
#colnames(pheno)[1] <- "scanID"

#gds <- GdsGenotypeReader(gdsfile)
#genoData <- GenotypeData(gds, scanAnnot=ScanAnnotationDataFrame(pheno))
#iterator <- GenotypeBlockIterator(genoData, snpBlock=variant_block_size, snpInclude=filtered)


test <- switch(tolower(test_type),
               score="Score",
               wald="Wald")


##updated function
assoc <- assocTestSingle(iterator, nullModel, test=test, imputed=TRUE)

##old assoc function
#assoc <- assocTestMM(seqData, nullModel,test=test)

#can add mac if desired
assoc <- addMAC(assoc, "single")

save(assoc, file=paste(out_prefix,".chr",chr,".RData",sep=""))

seqClose(seqData) #seqarray
#close(iterator)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
