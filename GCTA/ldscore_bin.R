argsv <- commandArgs(trailingOnly = TRUE)
chr <- argsv[1]

input.file <- paste0("chr", chr,".score.ld")

lds_seg = read.table(input.file,header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

write.table(lb1_snp, paste0("chr", chr,".snp_group1.txt"), row.names=F, quote=F, col.names=F)
write.table(lb2_snp, paste0("chr", chr,".snp_group2.txt"), row.names=F, quote=F, col.names=F)
write.table(lb3_snp, paste0("chr", chr,".snp_group3.txt"), row.names=F, quote=F, col.names=F)
write.table(lb4_snp, paste0("chr", chr,".snp_group4.txt"), row.names=F, quote=F, col.names=F)
