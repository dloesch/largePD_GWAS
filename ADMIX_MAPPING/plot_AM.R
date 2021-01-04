library(coda)

load("/mnt/hdd/large-PD/v2/AM/Admixture_mapping/Admixture_mapping_LARGE_PD_EUR.NAM_031220.RData")
chrs <- unique(EUR.NAM$chr[EUR.NAM$Joint.pval < 5E-04]) #just run on chromosomes with "interesting" results

#set significance level using empirical autoregression
tests <- c()

for(chr in 1:22){
  foo <- EUR.NAM[EUR.NAM$chr == chr,]
  fit <- spectrum0.ar(foo$Joint.pval)
  tests <- c(tests, fit$spec)
}
total <- sum(tests)
signif <- 0.05/total
print(signif)

#take whichever is more conservative
signif <- ifelse(signif > 5E-05, 5E-05, signif)

##### manhattan plot of AM results #####
data <- EUR.NAM

index <- 1:nrow(data)
data$index <- index
max_pval <- max(-log10(data$Joint.pval))
max_x <- ifelse(max_pval > signif, max_pval, signif)
max_x <- max_x

signif <- -log10(signif)


colors <- rep(c("slategray3","black"),times=11)
bins <- c()

chr <- 1
pval <- -log10(data$Joint.pval[data$chr == chr])
index <- data$index[data$chr == chr]
bins <- c(bins,median(index))


pdf("AM_resjuts.joint.pdf", width = 16, height = 8)
#png("AM_resjuts.joint.png",width = 1024, height = 768)
plot(pval~index, pch=20, xaxt="n", ylab="-log10 p-value", 
     xlab="chromosome", main="Admixture Mapping: Joint Test", type="p",
     ylim=c(0,max_x+0.5), xlim=c(0,nrow(data)+1), col=colors[chr])

for(chr in 2:22){
  pval <- -log10(data$Joint.pval[data$chr == chr])
  index <- data$index[data$chr == chr]
  bins <- c(bins,median(index))
  points(pval~index, pch=20, col=colors[chr])
}

#chromosome labels
axis(1, at=bins, labels=c(1:22), las=2)
#signif abline
abline(h=signif,lty=2, col="red")


dev.off()


####plot single ancestry results#####
load("Admixture_mapping/Admixture_mapping_LARGE_PD_AFR_031220.RData")
load("/mnt/hdd/large-PD/v2/AM/Admixture_mapping/Admixture_mapping_LARGE_PD_EUR_031220.RData")
load("Admixture_mapping/Admixture_mapping_LARGE_PD_NAM_031220.RData")

labels <- c("EUR", "NAM", "AFR")
AM <- list(EUR, NAM, AFR)

for(i in 1:3){
  data <- AM[[i]]
  
  index <- 1:nrow(data)
  data$index <- index
  max_pval <- max(-log10(data$pval))
  
  
  max_x <- ifelse(max_pval > signif, max_pval, signif)
  max_x <- max_x
  
  #manhattan plot
  colors <- rep(c("dodgerblue4","skyblue2"),times=11)
  bins <- c()
  
  chr <- 1
  pval <- -log10(data$pval[data$chr == chr])
  index <- data$index[data$chr == chr]
  bins <- c(bins,median(index))
  
  
  #png("AM_resjuts.joint.png",width = 1024, height = 768)
  plot(pval~index, pch=20, xaxt="n", ylab="-log10 p-value", 
       xlab="chromosome", main=paste("Admixture Mapping:", labels[i]), type="p",
       ylim=c(0,max_x+0.5), xlim=c(0,nrow(data)+1), col=colors[chr])
  
  for(chr in 2:22){
    pval <- -log10(data$pval[data$chr == chr])
    index <- data$index[data$chr == chr]
    bins <- c(bins,median(index))
    points(pval~index, pch=20, col=colors[chr])
  }
  
  #chromosome labels
  axis(1, at=bins, labels=c(1:22), las=2)
  #signif abline
  abline(h=signif,lty=2, col="firebrick1")
}

