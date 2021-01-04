library(SeqVarTools)
library(Biobase)
library(GENESIS)
library(SeqArray)

source("topmed.R")

#get file names
prefix <- "../GWAS/large"
null_model_file <- paste(prefix,".null_model.RData", sep="")
phenotype_file <- paste(prefix,".pheno.RData", sep="")


#gds file
gdsfile <- "../../replicate/large.Nalls_vars.gds"
gds <- seqOpen(gdsfile)

# createSeqVarData object
# get phenotypes
load(phenotype_file)
annot <- phen$annot
pData(annot)$sample.id <- as.character(pData(annot)$sample.id) #LRAGE-PD has numeric ids, change to character
#if genotype data has more subjects than phenotype data, use the next 3 lines
s <- as.data.frame(seqGetData(gds,"sample.id"))
colnames(s) <- "sample.id"
pData(annot) <- merge(pData(annot),s,by.x="sample.id",by.y="sample.id",all.x=TRUE,all.y=TRUE)
seqData <- SeqVarData(gds, sampleData=annot)

# get null model
nullModel <- getobj(null_model_file)

# get samples included in null model
sample.id <- nullModel$sample.id

checkSelectedVariants(seqData)

# create iterator, seqarray
block.size <- 1024
iterator <- SeqVarBlockIterator(seqData, variantBlock=block.size)

##updated function
assoc <- assocTestSingle(iterator, nullModel, test="Score", imputed=TRUE)

#add MAC
assoc <- addMAC(assoc, "single")

seqClose(seqData) #seqarray


###now compare with Nalls####
gds <- seqOpen(gds.fn = gdsfile)

#add SNP name, ref, alt
assoc$SNP <- seqGetData(gds, "annotation/id")
assoc$ref <- seqGetData(gds, "$ref")
assoc$alt <- seqGetData(gds, "$alt")

seqClose(gds)

#get betas
#convert score test to betas using this approximation to the Wald
assoc$beta <- assoc$Score/(assoc$Score.SE^2)

#get list of SNPS to update to rs numbers 
rsnum <- read.table("update_snps.txt", header=FALSE, stringsAsFactors = FALSE)
colnames(rsnum) <- c("SNP", "RSNUM")

#restrict to intersection of snps
l <- assoc[assoc$SNP %in% rsnum$SNP,]
l <- merge(l, rsnum, by.x="SNP", by.y="SNP", all.x=TRUE)

#read in Nalls et al. data
nalls <- read.table("PD.gwas_signif.NALLS.txt", header=TRUE, stringsAsFactors = FALSE)
l <- l[l$RSNUM %in% nalls$SNP,]

#set direction of effect base on effect alleles
for(i in 1:nrow(l)){
  if(l$alt[i] != nalls$EFFECT_ALLELE[nalls$SNP == l$RSNUM[i]]){
    l$BETA_MATCHED[i] <- l$beta[i]*-1
    l$AF_MATCHED[i] <- 1 - l$freq
  }else{
    l$BETA_MATCHED[i] <- l$beta[i]
    l$AF_MATCHED[i] <- l$freq[i]
  }
}

#prepare large-PD data
large <- l[c("RSNUM","chr", "pos", "ref", "alt", "BETA_MATCHED", "Score.pval", "AF_MATCHED", "MAC")]
colnames(large) <- c("SNP", "CHR", "POS","REF","ALT", "BETA_LARGE", "P", "AF_LARGE", "MAC")

#prepare nalls data
data <- nalls[c("SNP", "BETA")]
colnames(data)[2] <- "BETA_NALLS"

#prepapre output
data <- merge(data, large,by.x="SNP", by.y="SNP", all.x=FALSE)

#check for data concordance
data$BETA_CONCORDANCE <- ifelse(data$BETA_LARGE > 0 & data$BETA_NALLS >0, TRUE, 
  ifelse(data$BETA_LARGE < 0 & data$BETA_NALLS < 0, TRUE, FALSE))

#Bonferoni correction
data$P_ADJ <- data$P*nrow(data)
data$P_ADJ[data$P_ADJ > 1] <- 1

#classify by p-values
data$P_THRESH <- ifelse(data$P > 0.05, "Not significant", ifelse(data$P < 0.05, "Nominally Significant",NA))
data$P_THRESH[data$P < 0.05/nrow(data)] <- "Significant"
data$P_THRESH <- as.factor(data$P_THRESH)


#Add Nalls frequency data and nearest gene

nalls <- read.csv("./risk_variants.csv", header=TRUE, stringsAsFactors = FALSE)
nalls <- nalls[c(1,4,8)]
nalls <- nalls[complete.cases(af),]
colnames(nalls) <- c("SNP", "GENE", "AF_NALLS")

data <- merge(data, nalls, by.x="SNP", by.y="SNP", all.x=TRUE)

#set MAF from AF
data$MAF_NALLS <- ifelse(data$AF_NALLS > 0.5, 1-data$AF_NALLS, data$AF_NALLS)
data$MAF <- ifelse(data$AF > 0.5, 1- data$AF, data$AF)

#set filter for CG/AT sites, exclude if MAF > 0.3 and visually check allele frequencies to see if concordant
data$CGAT <- ifelse(data$REF == "G" & data$ALT == "C", TRUE,
  ifelse(data$REF == "C" & data$ALT == "G", TRUE, 
  ifelse(data$REF == "A" & data$ALT == "T", TRUE,
  ifelse(data$REF == "T" & data$ALT == "A", TRUE, FALSE))))

data$FILTER <- ifelse(data$MAF_NALLS > 0.3 & data$CGAT == TRUE, "FAIL", "PASS")
#no additional CG/AT sites need to be removed or flipped 

nrow(data[data$BETA_CONCORDANCE == TRUE & data$FILTER == "PASS",])/nrow(data[data$FILTER == "PASS",])


#set MAC filter, exclude if MAC <= 10
data$FILTER <- ifelse(data$MAC <=10, "FAIL", data$FILTER)

#find the N concordant sites
nrow(data[data$BETA_CONCORDANCE == TRUE & data$FILTER == "PASS",])
nrow(data[data$FILTER == "PASS",])
nrow(data[data$BETA_CONCORDANCE == TRUE & data$FILTER == "PASS",])/nrow(data[data$FILTER == "PASS",])
nrow(data[data$P_THRESH == "Nominally Significant",])

#subset data based on MAC and CG/AT filters
dat2 <- data[data$FILTER == "PASS",]
cor.test(dat2$BETA_NALLS, dat2$BETA_LARGE)

#set more restrictive MAF filter of 0.0452
dat3 <- dat2[dat2$MAF > 0.0452,]
cor.test(dat3$BETA_NALLS, dat3$BETA_LARGE)

#plot replication
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(ggrepel)


g <- ggscatter(data=dat2, x="BETA_LARGE", y="BETA_NALLS", color= "P_THRESH", add= "reg.line", add.params = list(color="black", alpha=0.4), 
                conf.int=FALSE, show.legend.text = TRUE)

g <- g +  geom_vline(xintercept = 0, color="darkgrey")
g <- g + geom_hline(yintercept = 0, color="darkgrey")

g <- g + annotate("rect", xmin=-Inf, xmax=0, ymin=-Inf, ymax=0, alpha=0.2, fill="slategray3")
g <- g + annotate("rect", xmin=0, xmax=Inf, ymin=0, ymax=Inf, alpha=0.2, fill="slategray3")


g <- g + labs(title = "Replication: CGAT sites removed", x = "LARGE-PD Betas",
                y="Nalls et al. Betas", color = "P-Value Thresholds")
g <- g+stat_cor(method = "pearson", label.x =0.2, label.y = 0.35) +theme_pubr()

g <- g + geom_label_repel(data = subset(dat2, P < 0.05), aes(label=as.character(GENE)),
                            box.padding   = 0.5, 
                            point.padding = 0.25,
                            label.padding =0.25,
                            segment.color = 'grey50',
                            segment.alpha = 0.6,
                            min.segment.length = 0.2,
                            direction="x",
                            nudge_y = -0.25,
                            size = 3)
g <- g+scale_color_manual(values=c("lightblue", "dodgerblue", "dodgerblue4"))


#plot 2, plot effect size difference (Nalls vs LARGE) by allele frequency
diff <- data
diff$diff <- abs(diff$BETA_NALLS - diff$BETA_LARGE)
diff$fill <- ifelse(diff$diff > (sd(diff$diff)+mean(diff$diff)),TRUE,FALSE)
diff$fill <- as.factor(diff$fill)
diff$AF <- ifelse(diff$AF > 0.5, 1-diff$AF, diff$AF)
d <- ggscatter(diff, x="AF", y="diff" , color="fill") +
  scale_color_brewer(palette = "Paired") +
  labs(title="Absolute beta differences by AF", color="> 1SD of mean beta diff", y= "abs( beta diff)")

d <- d + geom_vline(xintercept=max(diff$AF[diff$fill == TRUE]+0.004), linetype="dashed") +
  annotate("text", x=0.1, y=1.5, label=paste("MAF of",signif(max(diff$AF[diff$fill == TRUE]),digits = 2)))

#save plots
pdf("large-PD.GENESIS_replication.pdf", width=9)
print(g2)
print(d)
dev.off()

#sort by position and write out
data <- data[order(data$CHR, data$POS),]

write.table(data, "large-PD.GENESIS_replication.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names = TRUE)
