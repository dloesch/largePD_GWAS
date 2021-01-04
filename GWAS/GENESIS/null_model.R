library(Biobase)
library(GENESIS)
library(gdsfmt)
library(SeqVarTools)

#get TOPMed pipelne functions
source("/mnt/hdd/large-PD/v2/scripts/GENESIS/topmed.R")


#get prefix
argsv <- commandArgs(trailingOnly=T)
prefix <- argsv[1]

#GRM
grm_file <- argsv[2]
grm_type <- argsv[3]

# get the number of threads available
# this also sets MKL_NUM_THREADS, which should speed up matrix calculations if we are running parallel MKL
countThreads()

# get phenotypes
load(paste(prefix,".pheno.RData",sep=""))
annot <- phen[["annot"]]
outcome <- phen[["outcome"]]
covars <- phen[["covars"]]
sample.id <- phen[["sample.id"]]
group.var <- phen[["group.var"]]

#set family. Binomial or Gaussian
family <- binomial
#family <- gaussian

# kinship matrix or GRM
if(grm_type == "GCTA"){
  grm <- getGRM(grm_file, sample.id)
}else{
  load(paste(prefix,".pcrelate.RData",sep=""))
  grm <- pcrelateToMatrix(pcrel, scaleKin=2)
  grm_type <- "pcrelate"
}


## fit null model
hetvar <- FALSE
if(hetvar ==TRUE){
  ##allowing for heterogeneous variance between groups
  nullmod <- fitNullModel(annot, outcome=outcome, covars=covars,
                          cov.mat=grm, sample.id=sample.id,
                          family=family, group.var=group.var)
}else{
  #new function
  nullmod <- fitNullModel(annot, outcome=outcome, covars=covars,
                          cov.mat=grm, sample.id=sample.id,
                          family=family)
  #old function
  #nullmod <- fitNullMM(annot, outcome=outcome, covars=covars,
  #                      covMatList=grm, family=family)
}


#inverse normalization
inverse_normal <- FALSE
if(inverse_normal == TRUE){
  norm.option <- "by.group"
  rescale <- "residSD"
  
  nullmod <- fitNullInvNorm(nullmod, cov.mat=grm,
                            norm.option=norm.option, rescale=rescale)
}

if(grm_type == "GCTA"){
  save(nullmod, file=paste(prefix,".null_model.RData", sep=""))
}else{
  save(nullmod, file=paste(prefix,".null_model.pcrelate.RData", sep=""))
}


# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
