
library(SeqArray)
library(GMMAT)
library(ggplot2)

##get GWAS results

load("large.assoc_single.chr4.RData")
vars <- assoc$variant.id[assoc$Score.pval < 1E-07]

#identify GWAS significant snps to include as covariates
#subset gds file
gds.file <- "../../gds/large.chr4.gds"
gds <- seqOpen(gds.file)
seqSetFilter(gds, variant.id=vars)

#write subsetted file out to vcf file
seqGDS2VCF(gds,"chr4.GWAS_signif.vcf.gz")

#convert vcf to R-friendly format using plink
system("plink --vcf chr4.GWAS_signif.vcf.gz --recode A --out chr4.GWAS_signif")

#save gds of filtered file for future use
seqExport(gds, "chr4.GWAS_signif.gds")
seqClose(gds)
gds.file <- "chr4.GWAS_signif.gds"

#read in output of plink command
geno <- read.table("chr4.GWAS_signif.raw", header=TRUE)
geno <-geno[-c(1,3:6)]
colnames(geno) <- gsub("...._.$","",colnames(geno))

#get pheno and merge with geno
p <- read.csv("GWAS/large.pheno.03_2021.csv", header=TRUE)
p <- merge(p, geno, by.x="sample.id", by.y="IID", all.x=TRUE)
p$SEX <- as.factor(p$SEX)

#get GRM
GRM.file = "../GWAS/large.ALL.grm.RData"
load(GRM.file)

#convert to matrix and set rownames/colnames
GRM <- as.data.frame(grm$grm)
colnames(GRM) <- grm$sample.id
rownames(GRM) <- grm$sample.id
GRM <- as.matrix(GRM)

##first fit model without any snps
f=as.formula("PD_STATUS ~ AGE_MERGED + SEX + PC1 + PC2 + PC3 + PC4 + PC5")
null <- glmmkin(fixed=f, data=p, kins=GRM, id="sample.id")
glmm.score(null, infile=gds.file, outfile="no_conditions.txt")

test <- read.table("no_conditions.txt", header=TRUE, na.strings = ".")


#specify which snps to use as covariates
gwas <- read.delim("/mnt/hdd/large-PD/v2/GWAS/post_gwas/chr4.GWAS_signif.txt", stringsAsFactors = FALSE)
           paste0("chr4.",gwas$pos[gwas$SNP == rsnum[3]]))

# sart with rs356182
rsnum <- "rs356182"
conditions <- paste0("chr4",".",gwas$pos[gwas$SNP == rsnum])

###rs356182####
#fit null model for first conditional snp
f=as.formula(paste("PD_STATUS ~ AGE + SEX + PC1 + PC2 + PC3 + PC4 + PC5+",conditions[1]))

null1 <- glmmkin(fixed=f, data=p, kins=GRM, id="sample.id")

#now for score test
gds.file <- "chr4.GWAS_signif.gds"
glmm.score(null1, infile=gds.file, outfile="condition1.txt")

#read in score test results
test1 <- read.table("condition1.txt", header=TRUE, na.strings = ".")

#find snp with lowest p-value after adjusting for rs356182
conditions <- c(conditions,
                paste(unname(unlist(strsplit(as.character(test1$SNP[test1$PVAL == min(test1$PVAL)]), split = ":")))[1:2], collapse="."))
pos <- unname(unlist(strsplit(as.character(test1$SNP[test1$PVAL == min(test1$PVAL)]), split = ":")))[2]
rsnum <- c(rsnum, gwas$SNP[gwas$pos == pos])

snps <- as.character(test1$SNP[test1$PVAL < 0.05/nrow(test1)])

##adding snp 2 as a fixed effect##
f=as.formula(paste("PD_STATUS ~ AGE + SEX + PC1 + PC2 + PC3 + PC4 + PC5+",conditions[1], "+", conditions[2]))

null2 <- glmmkin(fixed=f, data=p, kins=GRM, id="sample.id")
glmm.score(null2, infile=gds.file, outfile="condition2.txt")

test2 <- read.table("condition2.txt", header=TRUE,na.strings = ".")




#compile results
results <- test
results[[paste0("P_",rsnum[1])]] <- test1$PVAL
results[[paste0("P_",rsnum[2])]]<- test2$PVAL
#results$P_rs356191 <- test3$PVAL

#adjust pvals
results$P_rs356182_adj <- ifelse(results$P_rs356182*nrow(results) < 1, results$P_rs356182*nrow(results),1)
results$P_rs356182_adjR <- ifelse(results$P_rs356182*223 < 1, results$P_rs356182*223, 1)
results$P_rs6830166_adj <- ifelse(results$P_rs6830166*nrow(results) < 1, results$P_rs6830166*nrow(results),1)
results$P_rs6830166_adjR <- ifelse(results$P_rs6830166*223 < 1, results$P_rs6830166*223,1)


#add rsnums to this file
ids <- gwas[c(1:3,6,7)]
ids$ID <- paste0("chr",ids$chr,":",ids$pos,":", ids$REF, ":", ids$ALT)
ids <- ids[c("SNP", "ID")]
colnames(ids)[1] <- "RSNUM"
results <- merge(results, ids, by.x="SNP", by.y="ID", all.x=TRUE)

#save file
write.table(results, "chr4.conditional_analysis.GMMAT.txt", 
            sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)


library(ggplot2)
pdf("chr4.conditional_analysis.pdf")
data <- results[complete.cases(results$P_rs356182),]
data$SNP <- as.character(data$SNP)
g <- ggplot(data=data, aes(1:nrow(data), y=-log10(data$P_rs356182)))+geom_point() +
  theme_classic()
g <- g + scale_x_discrete(name ="SNP", limits=c(data$SNP)) +
  theme(axis.text.x = element_text(size=8, angle=45, hjust=1))
g <- g + geom_hline(yintercept = -log10(0.05/nrow(results)), linetype="dashed", colour="red")
g <- g+geom_hline(yintercept = -log10(0.05/223), linetype="dashed")
g <- g+ labs(title = "Conditioning on rs356182", colour= "significant", y="-log10(p-palues)")

g

data <- results[complete.cases(results$P_rs6830166),]
data$SNP <- as.character(data$SNP)
g <- ggplot(data=data, aes(1:nrow(data), y=-log10(data$P_rs6830166)))+geom_point() +
  theme_classic()
g <- g + scale_x_discrete(name ="SNP", limits=c(data$SNP)) +
  theme(axis.text.x = element_text(size=8, angle=45, hjust=1))
g <- g + geom_hline(yintercept = -log10(0.05/nrow(results)), linetype="dashed", colour="red")
g <- g+geom_hline(yintercept = -log10(0.05/223), linetype="dashed")
g <- g+ labs(title = "Conditioning on rs356182 and rs6830166", colour= "significant", y="-log10(p-palues)")

g


dev.off()




#clean-up
system("rm condition1.txt condition2.txt condition3.txt no_conditions.txt")
