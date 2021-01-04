#!/usr/bin/bash
#$ -q threaded.q
#$ -pe thread 4
#$ -P toconnor-lab
#$ -l mem_free=10G
#$ -e LD.log
#$ -o LD.log
#$ -cwd
#$ -N LD.GCTA


chr=$1

#path to plink2
plink2=/usr/local/packages/plink-2.00.alpha/plink2

#path to gcta
gcta64=/usr/local/bin/gcta64

#path to genotype file
geno=/local/chib/toconnor_grp/LARGE-PD/imputed/chr$chr.dose.vcf.gz

#convert to plink
$plink2 --vcf $geno --maf 0.01 --out chr$chr --make-bed --memory 10000M

#calculate ld scores
$gcta64 --bfile chr$chr --ld-score-region 200 --out chr$chr --thread-num 4

#create bins by LD score:
/usr/local/bin/Rscript ldscore_bin.R $chr

#estimate GRM stratified by LD scre 
$gcta64 --bfile chr$chr --extract chr$chr.snp_group1.txt --make-grm --out chr$chr.test_group1 --thread-num 4
$gcta64 --bfile chr$chr --extract chr$chr.snp_group2.txt --make-grm --out chr$chr.test_group2 --thread-num 4
$gcta64 --bfile chr$chr --extract chr$chr.snp_group3.txt --make-grm --out chr$chr.test_group3 --thread-num 4
$gcta64 --bfile chr$chr --extract chr$chr.snp_group4.txt --make-grm --out chr$chr.test_group4 --thread-num 4
