library(Biobase)
library(tools)
source("/mnt/hdd/large-PD/v2/scripts/GENESIS/topmed.R")

#get original pheno file
argsv <- commandArgs(trailingOnly = TRUE)
pheno_file <- argsv[1]

#determine if csv or txt/tsv file
isCSV <- file_ext(pheno_file) == "csv"
if(isCSV){
  pheno <- read.csv(pheno_file, stringsAsFactors = FALSE, header=T, comment.char = "#")
}else{
  pheno <- read.table(pheno_file, stringsAsFactors = FALSE, header=T, comment.char = "#")
}

#get prefix 
prefix <- argsv[2]

#get covars
covars <- argsv[3]
foo <- strsplit(covars, split=",")
covars <- c()
for(i in 1:length(foo[[1]])){
  covars <- c(covars,foo[[1]][i])
}

#get number of pcs
n_pcs <- argsv[4]


#get outcome 

outcome <- argsv[5]


#create annotated data frame
#need sample id column to be named correctly.
colnames(pheno)[1] <- "sample.id"

#make sure sex is a factor
pheno$SEX <- as.factor(pheno$SEX)


#create dummary metadata, replace with actual metadata if desired. see example below
n <- ncol(pheno)
labels <- rep(NA,ncol(pheno))
metadata <- data.frame(labelDescriptions=labels)

#example metadata
#metadata <- data.frame(labelDescriptions=c("Subject's idetnifier", "Parkinson Status", "Subject's Age at Data Collection",
#              "Subject's Sex"))

annot <- AnnotatedDataFrame(pheno, metadata)


##Add PCA
pca_file <- paste(prefix,".pca.RData", sep="")
if (n_pcs > 0) {
  pca <- getobj(pca_file)
  pcs <- pca$vectors[,1:n_pcs,drop=FALSE]
  pccols <- paste0("PC", 1:n_pcs)
  colnames(pcs) <- pccols
  pcs <- pcs[match(annot$sample.id, rownames(pcs)),,drop=FALSE]
  rownames(pcs) <- NULL
  pData(annot) <- cbind(pData(annot), pcs)
} else {
  pccols <- NULL
}

covars <- c(covars, pccols)

#change if you have a category for allowing heterogenous variance
group.var <- NULL


#sample ids
sample.id <- pheno[,1]

#create list of phenotype elements
phen <- list(annot=annot, outcome=outcome, covars=covars, group.var=group.var, sample.id=sample.id)

outfile <- paste(prefix,".pheno.RData", sep="")
save(phen, file=outfile)
